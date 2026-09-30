// Recorded-frame timing, matching overlays, and cancellation of stale image loads.
const assert = require('node:assert/strict');
const {MapAnimationPlayer} = require('../../scripts/benchmarks/map_animation.js');
async function main() {
  const draws = [], timers = [], cancelled = [], frames = [];
  const canvas = {getContext: () => ({clearRect() {draws.length=0;}, drawImage(image) {draws.push(image.src);}})};
  const player = new MapAnimationPlayer(canvas, (frame,index,playing) => frames.push({frame,index,playing}), {
    load: async src => ({src,width:16,height:8}),
    schedule: (callback,delay) => {timers.push({callback,delay});return timers.length;},
    cancel: id => cancelled.push(id)
  });
  const states = [0,1,2].map(i => ({base:'map'+i,coast:'coast'+i,boundary:'boundary'+i,duration_ms:600,time:i*5}));
  await player.setFrames(states);
  assert.deepEqual(draws,['map2']);
  await player.setLayers(['boundary','coast']);
  assert.deepEqual(draws,['map2','boundary2','coast2']);
  await player.play();
  assert.deepEqual(draws,['map0','boundary0','coast0']);
  assert.equal(timers.at(-1).delay,600);
  await player.setSpeed(2);
  assert.equal(timers.at(-1).delay,300);
  timers.at(-1).callback(); await new Promise(resolve=>setImmediate(resolve));
  assert.deepEqual(draws,['map1','boundary1','coast1']);
  await player.stop();
  assert.deepEqual(draws,['map2','boundary2','coast2']);
  assert.equal(frames.at(-1).playing,false);
  assert.ok(cancelled.length>0);
  assert.throws(()=>player.setSpeed(0),/Invalid/);
  let resolveOld;
  player.load = src => src==='old'?new Promise(resolve=>{resolveOld=resolve;}):Promise.resolve({src,width:16,height:8});
  const old = player.setFrames([{base:'old'}]);
  await player.setFrames([{base:'new'}]);
  resolveOld({src:'old',width:16,height:8});await old;
  assert.deepEqual(draws,['new']);
  console.log('Map animation: synchronized overlays, speed, stop, and stale-load cancellation passed');
}
main().catch(error=>{console.error(error);process.exitCode=1;});
