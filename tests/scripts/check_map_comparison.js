const assert=require('node:assert/strict');
const {MapFieldComparison}=require('../../scripts/benchmarks/map_comparison.js');
async function main() {
  let requests=0;
  const arrays={a:[1,2,3,4],b:[1,-2,5,4],bad:[1,NaN,3,4]};
  const reader=new MapFieldComparison(async url=>{requests++;return new Float64Array(arrays[url]).buffer;});
  const a={url:'a',cells:4},b={url:'b',cells:4};
  const diff=await reader.compare(a,b);
  assert.deepEqual(diff,{changed_cells:2,mean_delta:-.5,mae:1.5,rmse:Math.sqrt(5),max_abs:4});
  assert.deepEqual(await reader.compare(b,a),{...diff,mean_delta:.5});
  assert.equal(requests,2);
  assert.deepEqual(await reader.compare(a,b,true),{changed_cells:2});
  await assert.rejects(reader.compare(a,{...b,cells:5}),/grid mismatch/);
  await assert.rejects(reader.field({url:'bad',cells:4}),/Non-finite/);
  await assert.rejects(reader.field({url:'bad',cells:5}),/size mismatch/);
  console.log('Map comparison: numeric precision, categorical labels, caching and invalid assets passed');
}
main().catch(error=>{console.error(error);process.exitCode=1;});
