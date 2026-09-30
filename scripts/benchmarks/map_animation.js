// Origin-agnostic playback of recorded raster frames and matching overlay layers.
// The compositor swaps complete frames atomically; it never interpolates fields.
class MapAnimationPlayer {
  constructor(canvas, onFrame, dependencies = {}) {
    this.canvas = canvas;
    this.onFrame = onFrame;
    this.load = dependencies.load || (src => new Promise((resolve, reject) => {
      const image = new Image();
      image.onload = () => resolve(image); image.onerror = () => reject(new Error('Cannot load map frame'));
      image.src = src;
    }));
    this.schedule = dependencies.schedule || ((callback, delay) => setTimeout(callback, delay));
    this.cancel = dependencies.cancel || (timer => clearTimeout(timer));
    this.cache = new Map(); this.frames = []; this.layers = []; this.speed = 1;
    this.index = 0; this.playing = false; this.timer = null; this.generation = 0;
  }
  image(src) {
    if (!this.cache.has(src)) this.cache.set(src, this.load(src));
    return this.cache.get(src);
  }
  async render() {
    if (this.timer !== null) this.cancel(this.timer);
    this.timer = null;
    const generation = ++this.generation, frame = this.frames[this.index];
    if (!frame) return;
    const sources = [frame.base, ...this.layers.map(layer => frame[layer]).filter(Boolean)];
    try {
      const images = await Promise.all(sources.map(src => this.image(src)));
      if (generation !== this.generation) return;
      this.canvas.width = images[0].width; this.canvas.height = images[0].height;
      const context = this.canvas.getContext('2d'); context.imageSmoothingEnabled = false;
      context.clearRect(0, 0, this.canvas.width, this.canvas.height);
      for (const image of images) context.drawImage(image, 0, 0);
      this.onFrame(frame, this.index, this.playing);
      if (this.playing) this.timer = this.schedule(() => {
        this.index = (this.index + 1) % this.frames.length; this.render();
      }, frame.duration_ms / this.speed);
    } catch (error) {
      if (generation !== this.generation) return;
      this.playing = false; this.onFrame(frame, this.index, false, error);
    }
  }
  setFrames(frames) {
    this.playing = false; this.frames = frames; this.index = Math.max(0, frames.length - 1);
    this.cache.clear(); return this.render();
  }
  setLayers(layers) { this.layers = layers; return this.render(); }
  setSpeed(speed) {
    if (!Number.isFinite(speed) || speed <= 0) throw new Error('Invalid animation speed');
    this.speed = speed; return this.render();
  }
  play() { this.playing = true; this.index = 0; return this.render(); }
  stop() { this.playing = false; this.index = Math.max(0, this.frames.length - 1); return this.render(); }
}
if (typeof module !== 'undefined') module.exports = {MapAnimationPlayer};
