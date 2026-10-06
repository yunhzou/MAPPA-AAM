/* Deterministic offline capture. PLAYWRIGHT_MODULE / CHROME / FFMPEG are overridable. */
const fs=require('fs'),path=require('path'),{spawn}=require('child_process');
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
async function main(){
 const dir=path.resolve(process.argv[2]||'output'),mode=process.argv[3]||'preview';
 const browser=await chromium.launch({executablePath:process.env.CHROME||'/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',headless:true,args:['--enable-unsafe-swiftshader','--disable-web-security']});
 const page=await browser.newPage({viewport:{width:1440,height:1024},deviceScaleFactor:1});
 const errors=[];page.on('pageerror',e=>errors.push(String(e)));
 await page.goto('file://'+dir+'/index.html');await page.waitForFunction(()=>!!window.film);await page.evaluate(()=>document.fonts.ready);
 await page.evaluate(()=>{const d=window.film.data;if(d.decoded_candidates.length!==2||new Set(d.decoded_candidates.map(x=>x.event_class)).size!==2)throw Error('Final event cards must be canonical-unique');if(d.paths.length!==2||new Set(d.paths.map(p=>p.event_class)).size!==2||d.paths.some(p=>p.event_class!==d.decoded_candidates[p.class_index].event_class))throw Error('Exactly two grown branches must decode to their distinct classes');});
 await page.evaluate(()=>{window.film.seek(19);const tree=document.querySelector('#tree');const labels=[...tree.querySelectorAll('[data-witness]')].map(n=>n.dataset.witness);if(labels.join(',')!=='A,B'||/B \+ C|Other saved path|Three displayed/.test(document.body.textContent))throw Error('Ending must show only grown witnesses A and B');});
 const duration=await page.evaluate(()=>window.film.duration);
 await page.evaluate(()=>{const a=window.film.seek(7.2),b=window.film.seek(8.2),c=window.film.seek(9.2);if(a.branchFocus.index!==0||b.branchFocus.index!==1||a.branchFocus.color===b.branchFocus.color||a.cameraTime!==b.cameraTime||c.path!==1)throw Error('Branch focus must switch color once while camera stays fixed');});
 const shot=async(t,file)=>{await page.evaluate(t=>window.film.seek(t),t);await page.locator('#film').screenshot({path:file,animations:'disabled'});};
 for(const t of [0,3,6,7.2,7.6,8.2,8.6,9.5,11.2,13.5,19])await shot(t,path.join(dir,`review-${t}.png`));
 // Verify every displayed frame can be selected and every mapped atom remains injective.
 for(let i=0;i<=240;i++){const t=i/240*duration;const state=await page.evaluate(t=>window.film.seek(t),t);if(new Set(Object.values(state.mapping)).size!==Object.keys(state.mapping).length)throw Error('Non-injective displayed mapping at '+t);}
 if(mode==='video'){
  const fps=24,encoder=spawn(process.env.FFMPEG||'ffmpeg',['-y','-f','image2pipe','-vcodec','png','-framerate',String(fps),'-i','-','-an','-c:v','libx264','-threads',process.env.FFMPEG_THREADS||'4','-preset','fast','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart',path.join(dir,'mappa-grow-branch-decode.mp4')],{stdio:['pipe','ignore','pipe']});
  let log='';encoder.stderr.on('data',b=>log+=b);const done=new Promise((resolve,reject)=>{encoder.on('close',code=>code?reject(Error(log)):resolve());encoder.on('error',reject)});
  for(let i=0;i<duration*fps;i++){
   await page.evaluate(t=>window.film.seek(t),i/fps);
   const buf=await page.locator('#film').screenshot({animations:'disabled'});
   if(!encoder.stdin.write(buf))await new Promise(resolve=>encoder.stdin.once('drain',resolve));
   if(i%(fps*4)===0)console.log(`${i/fps}/${duration} seconds rendered`);
  }encoder.stdin.end();await done;
 }
 fs.writeFileSync(path.join(dir,'browser-validation.json'),JSON.stringify({status:errors.length?'failed':'passed',errors,preview_times:[0,3,6,7.2,7.6,8.2,8.6,9.5,11.2,13.5,19],branch_focus_check:true,final_witness_labels:["A","B"],final_event_class_count:2,sampled_timeline_injections:241,viewport:[1440,900],mode},null,2));
 await browser.close();if(errors.length)throw Error(errors.join('\n'));
}
main().catch(e=>{console.error(e);process.exit(1)});
