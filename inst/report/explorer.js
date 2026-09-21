/* SpaMTP offline region explorer. No remote libraries or data requests. */
(() => {
  'use strict';
  const node = document.getElementById('spamtp-explorer-data');
  if (!node) return;
  const data = JSON.parse(new TextDecoder().decode(Uint8Array.from(atob(JSON.parse(node.textContent)), c => c.charCodeAt(0))));
  const el = id => document.getElementById('ex-' + id);
  const finite = x => typeof x === 'number' && Number.isFinite(x);
  const fmt = x => finite(x) ? (Math.abs(x) > 1e4 || (x !== 0 && Math.abs(x) < .001) ? x.toExponential(3) : Number(x.toPrecision(4)).toString()) : '—';
  const state = {modality: 0, mode: 0, feature: null, page: 0, sort: 'default', ascending: false, zoom: 1, pan: [0, 0], filtered: []};
  let effectHits = [], spatialHits = [], distributionHits = [];
  const tooltip = el('tooltip');
  const mod = () => data.modalities[state.modality];
  const mode = () => mod().modes[state.mode] || null;
  const profile = () => mod().profiles.find(p => p.feature === state.feature);
  function options(select, choices) {
    select.replaceChildren();
    choices.forEach(([value, label]) => { const option = document.createElement('option'); option.value = value; option.textContent = label; select.append(option); });
  }
  function context(id) {
    const canvas = el(id), ctx = canvas.getContext('2d');
    ctx.clearRect(0, 0, canvas.width, canvas.height); ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.font = '12px system-ui'; ctx.strokeStyle = '#d7e1e7'; ctx.fillStyle = '#304b5b';
    return {canvas, ctx, width: canvas.width, height: canvas.height};
  }
  function extent(values) {
    const numbers = values.filter(finite); if (!numbers.length) return [-1, 1];
    let lo = Infinity, hi = -Infinity; numbers.forEach(v => {lo = Math.min(lo, v); hi = Math.max(hi, v);});
    if (lo === hi) {lo -= .5; hi += .5;} return [lo, hi];
  }
  function axes(c, xlim, ylim, xlabel, ylabel) {
    const {ctx, width, height} = c, left = 82, top = 24, right = width - 20, bottom = height - 52;
    const x = v => left + (v - xlim[0]) / (xlim[1] - xlim[0]) * (right - left);
    const y = v => bottom - (v - ylim[0]) / (ylim[1] - ylim[0]) * (bottom - top);
    ctx.lineWidth = 1; ctx.strokeStyle = '#aebfc9'; ctx.beginPath(); ctx.moveTo(left, top); ctx.lineTo(left, bottom); ctx.lineTo(right, bottom); ctx.stroke();
    ctx.fillStyle = '#526977'; ctx.textAlign = 'center';
    for (let i = 0; i <= 4; i++) {let v = xlim[0] + (xlim[1]-xlim[0])*i/4; ctx.fillText(fmt(v), x(v), bottom + 20);}
    ctx.textAlign = 'right'; for (let i = 0; i <= 4; i++) {let v = ylim[0]+(ylim[1]-ylim[0])*i/4;ctx.fillText(fmt(v),left-8,y(v)+4);}
    ctx.textAlign = 'center'; ctx.fillText(xlabel,width/2,height-10);
    ctx.save();ctx.translate(14,height/2);ctx.rotate(-Math.PI/2);ctx.fillText(ylabel,0,0);ctx.restore();
    return {x, y, left, right, top, bottom};
  }
  function dot(ctx, x, y, colour, radius=3) {ctx.fillStyle=colour;ctx.beginPath();ctx.arc(x,y,radius,0,Math.PI*2);ctx.fill();}
  function message(c, text) {c.ctx.fillStyle='#607986';c.ctx.textAlign='center';c.ctx.fillText(text,c.width/2,c.height/2);}
  function showTip(event, lines) {tooltip.textContent=lines.join('\n');tooltip.hidden=false;tooltip.style.left=Math.min(event.clientX+12,window.innerWidth-330)+'px';tooltip.style.top=Math.min(event.clientY+12,window.innerHeight-120)+'px';}
  function mousePoint(event, canvas) {const rect=canvas.getBoundingClientRect();return [(event.clientX-rect.left)*canvas.width/rect.width,(event.clientY-rect.top)*canvas.height/rect.height];}
  function nearest(event, canvas, hits) {const p=mousePoint(event,canvas);let best=null,d=81;hits.forEach(h=>{const v=(h.x-p[0])**2+(h.y-p[1])**2;if(v<d){d=v;best=h;}});return best;}
  function filter() {
    const m=mode();if(!m){state.filtered=[];return;}
    const search=el('search').value.toLowerCase(), direction=el('direction').value;
    const threshold=Math.max(0,Number(el('effect').value)||0), cutoff=Math.max(0,Math.min(1,Number(el('fdr').value)));
    state.filtered=m.rows.filter(r=>finite(r[1]) && Math.abs(r[1])>=threshold &&
      (direction==='all'||(direction==='up'?r[1]>0:r[1]<0)) &&
      (!search||mod().labels[r[0]].toLowerCase().includes(search)) &&
      (m.kind!=='test'||(finite(r[6])&&r[6]<=cutoff)) &&
      (m.kind==='test'||!finite(r[8])||(direction==='up'?r[8]:direction==='down'?1-r[8]:Math.max(r[8],1-r[8]))>=Number(el('auc').value)) &&
      (m.kind==='test'||Number(el('stability').value)===0||(finite(r[12])&&r[12]>=Number(el('stability').value))));
    if(state.sort==='feature')state.filtered.sort((a,b)=>mod().features[a[0]].localeCompare(mod().features[b[0]])*(state.ascending?1:-1));
    if(state.sort==='effect')state.filtered.sort((a,b)=>(a[1]-b[1])*(state.ascending?1:-1));
    if(state.sort==='auc')state.filtered.sort((a,b)=>((a[8]??0)-(b[8]??0))*(state.ascending?1:-1));
    if(state.sort==='cohen')state.filtered.sort((a,b)=>((a[10]??0)-(b[10]??0))*(state.ascending?1:-1));
    if(state.sort==='fdr')state.filtered.sort((a,b)=>((a[6]??Infinity)-(b[6]??Infinity))*(state.ascending?1:-1));
    state.page=Math.min(state.page,Math.max(0,Math.ceil(state.filtered.length/20)-1));
  }
  function effects() {
    const c=context('effects'),m=mode();effectHits=[];
    if(!m){message(c,'No eligible region or replicate contrast.');return;}
    if(!state.filtered.length){message(c,'No features match these filters.');return;}
    const tested=m.kind==='test', ys=state.filtered.map(r=>tested?-Math.log10(Math.max(r[6]??1,1e-300)):finite(r[8])?r[8]:(r[4]??0)-(r[5]??0));
    const ylim=tested?[0,ys.reduce((v,y)=>Math.max(v,y),1)]:extent(ys), xlim=extent(state.filtered.map(r=>r[1]));
    const a=axes(c,xlim,ylim,'Effect (workflow units)',tested?'-log10(FDR)':state.filtered.some(r=>finite(r[8]))?'Mean pairwise AUC':'Non-zero fraction difference');
    if(xlim[0]<0&&xlim[1]>0){c.ctx.strokeStyle='#d3dde2';c.ctx.beginPath();c.ctx.moveTo(a.x(0),a.top);c.ctx.lineTo(a.x(0),a.bottom);c.ctx.stroke();}
    state.filtered.forEach((r,i)=>{const h={x:a.x(r[1]),y:a.y(ys[i]),row:r,feature:r[0]};effectHits.push(h);dot(c.ctx,h.x,h.y,r[1]>=0?'#b75442aa':'#337da4aa',r[0]===state.feature?5:2.8);});
  }
  function spatial() {
    const c=context('spatial');spatialHits=[];
    if(!data.coordinates){message(c,'No spatial coordinates supplied for these observations.');return;}
    const p=profile(),m=mode(),chosen=el('sample').value;
    const indexes=data.coordinates.map((_,i)=>i).filter(i=>chosen==='all'||data.samples[i]===chosen);
    if(!indexes.length){message(c,'No preview observations in this sample.');return;}
    const xs=extent(indexes.map(i=>data.coordinates[i][0])),ys=extent(indexes.map(i=>data.coordinates[i][1]));
    const scale=Math.min((c.width-70)/(xs[1]-xs[0]),(c.height-60)/(ys[1]-ys[0]))*state.zoom;
    const cx=(xs[0]+xs[1])/2,cy=(ys[0]+ys[1])/2,range=p?extent(indexes.map(i=>p.values[i])):[0,1];
    c.ctx.save();c.ctx.beginPath();c.ctx.rect(5,5,c.width-10,c.height-40);c.ctx.clip();
    indexes.forEach(i=>{
      const x=(data.coordinates[i][0]-cx)*scale+c.width/2+state.pan[0],y=(data.coordinates[i][1]-cy)*scale+(c.height-35)/2+state.pan[1];
      let colour='#c8d3da';
      if(el('colour').value==='expression'&&p){const t=Math.max(0,Math.min(1,(p.values[i]-range[0])/(range[1]-range[0])));colour=`rgb(${Math.round(235-65*t)},${Math.round(240-185*t)},${Math.round(243-184*t)})`;}
      else if(m){const label=m.kind==='test'&&data.conditions?data.conditions[i]:data.regions[i];colour=label===m.target?'#ba5746':'#c8d3da';}
      dot(c.ctx,x,y,colour,2.4);spatialHits.push({x,y,index:i});
    });c.ctx.restore();
    c.ctx.fillStyle='#526977';c.ctx.textAlign='left';
    c.ctx.fillText(el('colour').value==='expression'&&p?`Expression: ${fmt(range[0])} to ${fmt(range[1])}`:'Selected region highlighted',15,c.height-12);
    c.ctx.textAlign='right';c.ctx.fillText(indexes.length+' preview observations',c.width-12,c.height-12);
  }
  function means() {
    const c=context('means'),d=mod();
    if(state.feature===null||!d.region_means||!d.region_names.length){message(c,'Region summaries unavailable.');return;}
    const all=d.region_names.map((n,i)=>({name:n,value:d.region_means[state.feature][i],index:i,count:d.region_counts?.[i]})).filter(r=>finite(r.value));
    all.sort((a,b)=>b.value-a.value);const selected=all.slice(0,12),lim=extent([0,...selected.map(r=>r.value)]);
    const ctx=c.ctx,left=150,right=c.width-85,step=(c.height-35)/Math.max(1,selected.length);
    const x=v=>left+(v-lim[0])/(lim[1]-lim[0])*(right-left);
    selected.forEach((r,i)=>{const y=12+i*step;ctx.fillStyle=mode()&&mode().target===r.name?'#ba5746':'#448a9e';ctx.fillRect(Math.min(x(0),x(r.value)),y,Math.max(1,Math.abs(x(r.value)-x(0))),Math.min(18,step-3));ctx.fillStyle='#3c5665';ctx.textAlign='right';ctx.fillText(r.name.slice(0,16)+(finite(r.count)?' (n='+r.count+')':''),left-8,y+12);ctx.textAlign='left';ctx.fillText(fmt(r.value),x(r.value)+4,y+12);});
    if(all.length>12){ctx.textAlign='left';ctx.fillText('12 highest regional means shown; complete values are in the CSV.',12,c.height-5);}
  }
  function distribution() {
    const c=context('distribution'),p=profile();distributionHits=[];
    if(!p){message(c,'Spatial/distribution profile is not embedded for this feature.');return;}
    const labels=[...new Set(data.regions.filter(g=>g!==null))],m=mode();
    const ordered=labels.slice().sort((a,b)=>(a===m?.target?-1:b===m?.target?1:0)).slice(0,8);
    const indexes=p.values.map((_,i)=>i).filter(i=>ordered.includes(data.regions[i])&&(el('sample').value==='all'||data.samples[i]===el('sample').value));
    const ylim=extent(indexes.map(i=>p.values[i]));const a=axes(c,[0,Math.max(1,ordered.length)],ylim,'Regions (up to 8 in preview)','Workflow value');
    c.ctx.fillStyle='white';c.ctx.fillRect(30,c.height-43,c.width-30,24);
    indexes.forEach(i=>{const j=ordered.indexOf(data.regions[i]),jitter=((i*2654435761>>>0)%1000)/1000;const x=a.x(j+.2+.6*jitter),y=a.y(p.values[i]);dot(c.ctx,x,y,data.regions[i]===m?.target?'#b7544280':'#3e819c50',1.8);distributionHits.push({x,y,index:i});});
    c.ctx.fillStyle='#526977';c.ctx.textAlign='center';ordered.forEach((name,i)=>c.ctx.fillText(String(name).slice(0,10),a.x(i+.5),c.height-29));
  }
  function table() {
    const body=el('table').querySelector('tbody');body.replaceChildren();
    const start=state.page*20;
    state.filtered.slice(start,start+20).forEach(r=>{
      const tr=document.createElement('tr');tr.classList.toggle('selected',r[0]===state.feature);
      const td=document.createElement('td'),button=document.createElement('button');button.type='button';button.className='feature-button';button.textContent=mod().labels[r[0]];button.title=mod().labels[r[0]];button.onclick=()=>selectFeature(r[0]);td.append(button);tr.append(td);
      [1,2,3,4,5,8,10,12,6].forEach(i=>{const t=document.createElement('td');t.textContent=fmt(r[i]);tr.append(t);});body.append(tr);
    });
    const m=mode();el('count').textContent=`${state.filtered.length} matching / ${m?.rows.length??0} preloaded / ${m?.total??0} total; page ${state.page+1}`;
    el('prev').disabled=state.page===0;el('next').disabled=start+20>=state.filtered.length;
  }
  function detail() {
    el('feature-title').textContent=state.feature===null?'Select a feature':mod().labels[state.feature];
    const p=profile();el('profile-note').textContent=`All labelled observations contribute to regional means. ${data.observations.length}/${data.total_observations} observations are embedded for spatial/distribution browsing; ${mod().profiles.length}/${mod().features.length} feature profiles are preloaded. `+(p?'Click points for values; wheel to zoom and drag the spatial map.':'This feature has regional means and DE statistics; increase preview_features when rendering to embed its spatial profile.');
    means();distribution();spatial();
  }
  function selectFeature(f) {state.feature=f;detail();effects();table();}
  function refresh() {
    filter();const m=mode();
    el('effect-title').textContent=m?.kind==='test'?'DE volcano (biological replicates)':'Region effects and marker specificity';
    el('fdr').disabled=m?.kind!=='test';
    el('auc').disabled=m?.kind==='test';el('stability').disabled=m?.kind==='test';
    el('note').textContent=m?`${m.method} Displaying a preloaded ranking of ${m.rows.length}/${m.total} features; filters apply to this ranking. Full tables are linked below.`:'No eligible region comparisons. Choose a region column or inspect the native feature maps below.';
    if(state.feature===null||!mod().features[state.feature])state.feature=state.filtered.find(r=>mod().profiles.some(p=>p.feature===r[0]))?.[0]??mod().profiles[0]?.feature??null;
    effects();table();detail();
  }
  function changeMode() {state.mode=Number(el('mode').value)||0;state.page=0;state.sort='default';state.feature=null;el('fdr').value='.05';refresh();}
  function changeModality() {state.modality=Number(el('modality').value)||0;options(el('mode'),mod().modes.map((m,i)=>[i,m.title]));state.mode=0;state.feature=null;state.page=0;el('search').value='';changeMode();}
  options(el('modality'),data.modalities.map((m,i)=>[i,m.name]));
  options(el('sample'),[['all','All samples'],...[...new Set(data.samples)].map(s=>[s,s])]);
  el('modality').onchange=changeModality;el('mode').onchange=changeMode;
  ['search','direction','effect','auc','stability','fdr'].forEach(id=>{el(id).addEventListener('input',()=>{state.page=0;refresh();});});
  el('sample').onchange=()=>{state.zoom=1;state.pan=[0,0];detail();};el('colour').onchange=spatial;
  el('reset').onclick=()=>{state.zoom=1;state.pan=[0,0];spatial();};
  el('prev').onclick=()=>{state.page--;table();};el('next').onclick=()=>{state.page++;table();};
  el('table').querySelectorAll('[data-sort]').forEach(button=>{button.onclick=()=>{state.ascending=state.sort===button.dataset.sort?!state.ascending:true;state.sort=button.dataset.sort;refresh();};});
  el('effects').onclick=event=>{const hit=nearest(event,el('effects'),effectHits);if(hit)selectFeature(hit.feature);};
  el('effects').onmousemove=event=>{const hit=nearest(event,el('effects'),effectHits);if(hit)showTip(event,[mod().labels[hit.feature],'Effect: '+fmt(hit.row[1]),'AUC: '+fmt(hit.row[8]),'Cohen effect: '+fmt(hit.row[10]),'Direction stability: '+fmt(hit.row[12]),'FDR: '+fmt(hit.row[6])]);else tooltip.hidden=true;};
  ['spatial','distribution'].forEach(id=>{el(id).onmousemove=event=>{const hit=nearest(event,el(id),id==='spatial'?spatialHits:distributionHits);if(hit)showTip(event,[data.observations[hit.index],`Region: ${data.regions[hit.index]} | Sample: ${data.samples[hit.index]}`,'Value: '+fmt(profile()?.values[hit.index])]);else tooltip.hidden=true;};el(id).onclick=event=>{
    if(dragged){dragged=false;return;}const hit=nearest(event,el(id),id==='spatial'?spatialHits:distributionHits);
    if(hit){if(id==='spatial'&&el('colour').value==='region'){const index=mod().modes.findIndex(m=>m.kind==='region'&&m.target===data.regions[hit.index]);if(index>=0){el('mode').value=String(index);changeMode();}}
      showTip(event,[data.observations[hit.index],String(data.regions[hit.index]),'Value: '+fmt(profile()?.values[hit.index])]);}
  };});
  ['effects','spatial','distribution'].forEach(id=>el(id).addEventListener('mouseleave',()=>tooltip.hidden=true));
  let drag=null,dragged=false;el('spatial').addEventListener('pointerdown',e=>{drag=mousePoint(e,el('spatial'));dragged=false;el('spatial').setPointerCapture(e.pointerId);});
  el('spatial').addEventListener('pointermove',e=>{if(!drag)return;const p=mousePoint(e,el('spatial'));if(Math.abs(p[0]-drag[0])+Math.abs(p[1]-drag[1])>3)dragged=true;state.pan[0]+=p[0]-drag[0];state.pan[1]+=p[1]-drag[1];drag=p;spatial();});
  el('spatial').addEventListener('pointerup',()=>{drag=null;});el('spatial').addEventListener('pointercancel',()=>{drag=null;});
  el('spatial').addEventListener('wheel',e=>{e.preventDefault();state.zoom=Math.max(.5,Math.min(20,state.zoom*(e.deltaY<0?1.15:1/1.15)));spatial();},{passive:false});
  el('download').onclick=()=>{
    const quote=x=>'"'+String(x??'').replaceAll('"','""')+'"';
    const output=[['feature','effect','mean_region','mean_reference','detected_region','detected_reference','FDR','P.Value','mean_auc','min_auc','mean_cohen','min_cohen','direction_stability','rank_worst_block','cohen_min_block','cohen_max_block'],...state.filtered.map(r=>[mod().features[r[0]],...r.slice(1)])].map(r=>r.map(quote).join(',')).join('\n');
    const url=URL.createObjectURL(new Blob([output],{type:'text/csv;charset=utf-8'})),a=document.createElement('a');a.href=url;a.download='spamtp_filtered_features.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  };
  window.spaMTPExplorer={getState:()=>({...state,filteredCount:state.filtered.length,profileAvailable:!!profile()}),selectFeature, data};
  changeModality();document.getElementById('spamtp-explorer').dataset.ready='true';
})();
