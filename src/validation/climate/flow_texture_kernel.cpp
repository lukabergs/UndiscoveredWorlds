// Headless rendering accelerator. Numerical climate fields are never modified.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#ifdef _WIN32
#define UW_TEXTURE_API extern "C" __declspec(dllexport)
#else
#define UW_TEXTURE_API extern "C"
#endif
namespace {
constexpr double pi = 3.14159265358979323846;
constexpr double radius = 6371000.0;
double wrap(double x, int w) { return x - std::floor(x / w) * w; }
double sample(const double* a, int w, int h, double x, double y) {
    x = wrap(x, w); y = std::clamp(y, 0.0, double(h - 1));
    const int x0 = int(std::floor(x)), y0 = int(std::floor(y));
    const int x1 = (x0 + 1) % w, y1 = std::min(y0 + 1, h - 1);
    const double fx = x - x0, fy = y - y0;
    const double weights[4] = {(1-fx)*(1-fy), fx*(1-fy), (1-fx)*fy, fx*fy};
    const int indices[4] = {y0*w+x0, y0*w+x1, y1*w+x0, y1*w+x1};
    double result = 0;
    for (int i = 0; i < 4; ++i) if (weights[i] > 0) result += a[indices[i]] * weights[i];
    return result;
}
bool velocity(const double* east, const double* north, int w, int h, double x, double y, double& dx, double& dy) {
    const double latitude = (90.0 - (y + .5) * 180.0 / h) * pi / 180.0;
    dx = sample(east,w,h,x,y) / (2*pi*radius/w*std::max(std::cos(latitude),1e-6));
    dy = -sample(north,w,h,x,y) / (pi*radius/h);
    return std::isfinite(dx) && std::isfinite(dy);
}
bool direction(const double* east, const double* north, int w, int h, double x, double y, double& dx, double& dy) {
    if (!velocity(east,north,w,h,x,y,dx,dy)) return false;
    const double length = std::hypot(dx,dy);
    if (length <= 1e-12) return false;
    dx /= length; dy /= length; return true;
}
}
UW_TEXTURE_API void uw_wind_lic(const double* east, const double* north, const double* noise, int w, int h, float* output) {
    for (int row = 0; row < h; ++row) for (int col = 0; col < w; ++col) {
        const int cell = row*w+col;
        if (!std::isfinite(east[cell]) || !std::isfinite(north[cell])) {
            output[cell] = std::numeric_limits<float>::quiet_NaN(); continue;
        }
        double total = noise[cell], weightSum = 1;
        for (int sign : {-1,1}) {
            double x = col, y = row;
            for (int step = 1; step <= 12; ++step) {
                double dx,dy;
                if (!direction(east,north,w,h,x,y,dx,dy)) break;
                const double mx = wrap(x+sign*dx*.325,w), my = y+sign*dy*.325;
                if (my < 0 || my > h-1 || !direction(east,north,w,h,mx,my,dx,dy)) break;
                x = wrap(x+sign*dx*.65,w); y += sign*dy*.65;
                if (y < 0 || y > h-1) break;
                const double weight = .5+.5*std::cos(pi*step/13);
                total += weight*sample(noise,w,h,x,y); weightSum += weight;
            }
        }
        output[cell] = float(std::clamp((total/weightSum-.5)*2.2+.5,0.0,1.0));
    }
}
UW_TEXTURE_API void uw_wind_particles(const double* east, const double* north, int w, int h,
    const double* seedX, const double* seedY, int count, double seconds, double* hits) {
    std::fill(hits,hits+std::size_t(w)*h,0.0);
    for (int particle = 0; particle < count; ++particle) {
        double x = seedX[particle], y = seedY[particle], remaining = seconds;
        while (remaining > 0) {
            double dx,dy;
            if (!velocity(east,north,w,h,x,y,dx,dy)) break;
            const double rate = std::hypot(dx,dy);
            if (rate <= 1e-12) break;
            const double dt = std::min(remaining,std::min(10800.0,.5/rate));
            const double mx = wrap(x+dx*dt/2,w), my = y+dy*dt/2;
            if (my < 0 || my > h-1 || !velocity(east,north,w,h,mx,my,dx,dy)) break;
            const double nx = wrap(x+dx*dt,w), ny = y+dy*dt;
            if (!std::isfinite(nx) || !std::isfinite(ny) || ny < 0 || ny > h-1) break;
            x = nx; y = ny; remaining -= dt;
            hits[int(std::nearbyint(y))*w+(int(std::nearbyint(x))%w)] += 1;
        }
    }
}
