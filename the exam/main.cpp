#include <bits/stdc++.h>
using namespace std;

struct Point {
    double x, y, yp;
};

struct CubicHermiteInterval {
    double x0, x1, y0, y1, yp0, yp1;
    double h;

    static double H00(double t){ return 2*t*t*t - 3*t*t + 1; }
    static double H01(double t){ return -2*t*t*t + 3*t*t; }
    static double H10(double t){ return t*t*t - 2*t*t + t; }
    static double H11(double t){ return t*t*t - t*t; }

    double eval(double x) const {
        double t = (x - x0) / h;
        return y0*H00(t) + y1*H01(t) + (h*yp0)*H10(t) + (h*yp1)*H11(t);
    }

    double evalDer(double x) const {
        double t = (x - x0) / h;

        auto dH00 = [](double t){ return 6*t*t - 6*t; };
        auto dH01 = [](double t){ return -6*t*t + 6*t; };
        auto dH10 = [](double t){ return 3*t*t - 4*t + 1; };
        auto dH11 = [](double t){ return 3*t*t - 2*t; };

        double dtdx = 1.0 / h;
        double d_by_dt = y0*dH00(t) + y1*dH01(t) + (h*yp0)*dH10(t) + (h*yp1)*dH11(t);
        return d_by_dt * dtdx;
    }

    double evalSecondDer(double x) const {
        double t = (x - x0) / h;

        auto d2H00 = [](double t){ return 12*t - 6; };
        auto d2H01 = [](double t){ return -12*t + 6; };
        auto d2H10 = [](double t){ return 6*t - 4; };
        auto d2H11 = [](double t){ return 6*t - 2; };

        double d2_by_dt2 = y0*d2H00(t) + y1*d2H01(t) + (h*yp0)*d2H10(t) + (h*yp1)*d2H11(t);
        return d2_by_dt2 / (h*h);
    }
};

struct CubicSpline {
    vector<Point> pts;
    vector<CubicHermiteInterval> segs;

    void build(const vector<Point>& p) {
        pts = p;
        int n = (int)pts.size();
        segs.resize(n-1);
        for(int i=0;i<n-1;i++){
            segs[i] = {pts[i].x, pts[i+1].x, pts[i].y, pts[i+1].y, pts[i].yp, pts[i+1].yp, pts[i+1].x-pts[i].x};
        }
    }

    int findSeg(double x) const {
        int n = (int)pts.size();
        if(x <= pts[0].x) return 0;
        if(x >= pts[n-1].x) return n-2;
        int lo=0, hi=n-2;
        while(lo<=hi){
            int mid=(lo+hi)/2;
            if(x < pts[mid].x) hi=mid-1;
            else if(x > pts[mid+1].x) lo=mid+1;
            else return mid;
        }
        return max(0,min(n-2,lo));
    }

    double eval(double x) const { return segs[findSeg(x)].eval(x); }
    double evalDer(double x) const { return segs[findSeg(x)].evalDer(x); }
    double evalSecondDer(double x) const { return segs[findSeg(x)].evalSecondDer(x); }
};

static inline double quarticBasis(double x, double x0, double x1) {
    // (x-x0)^2 (x1-x)^2
    double a = x - x0;
    double b = x1 - x;
    return a*a*b*b;
}

struct QuarticSpline {
    CubicSpline base; // uses the cubic Hermite
    vector<double> e; // one per interval

    void build(const vector<Point>& p, bool e0_is_zero=true) {
        base.build(p);
        int n = (int)p.size();
        int m = n-1;
        e.assign(m, 0.0);

        if(!e0_is_zero) {
            // default is e[0]=0; adjust if you want different closure
        }

        // Enforce S'' continuity at interior knots x_i for i=1..n-2:
        // base left C'' + e[left]*Q''left == base right C'' + e[right]*Q''right
        // Solve for e[right] = e[i]
        for(int i=1;i<=n-2;i++){
            double xk = p[i].x;
            int leftInterval = i-1;  // [x_{i-1}, x_i]
            int rightInterval = i;  // [x_i, x_{i+1}]

            double C_left = base.segs[leftInterval].evalSecondDer(xk);
            double C_right = base.segs[rightInterval].evalSecondDer(xk);

            double hL = p[i].x - p[i-1].x;
            double hR = p[i+1].x - p[i].x;

            double Qpp_left = 2.0*hL*hL;
            double Qpp_right = 2.0*hR*hR;

            double rhs = C_left + e[leftInterval]*Qpp_left - C_right;
            e[rightInterval] = rhs / Qpp_right;
        }
    }

    int findSeg(double x) const {
        int n = (int)base.pts.size();
        if(x <= base.pts[0].x) return 0;
        if(x >= base.pts[n-1].x) return n-2;
        int lo=0, hi=n-2;
        while(lo<=hi){
            int mid=(lo+hi)/2;
            if(x < base.pts[mid].x) hi=mid-1;
            else if(x > base.pts[mid+1].x) lo=mid+1;
            else return mid;
        }
        return max(0,min(n-2,lo));
    }

    double eval(double x) const {
        int i = findSeg(x);
        auto &s = base.segs[i];
        double q = quarticBasis(x, s.x0, s.x1);
        return s.eval(x) + e[i]*q;
    }

    double evalDer(double x) const {
        int i = findSeg(x);
        auto &s = base.segs[i];
        double x0=s.x0, x1=s.x1;
        double a = x - x0;
        double b = x1 - x;

        // q'(x) for q=a^2 b^2
        double qprime = 2*a*b*b - 2*a*a*b;
        return s.evalDer(x) + e[i]*qprime;
    }

    double evalSecondDer(double x) const {
        int i = findSeg(x);
        auto &s = base.segs[i];
        double x0=s.x0, x1=s.x1;
        double a = x - x0;
        double b = x1 - x;

        // q'' for q=a^2 b^2 (with da/dx=1, db/dx=-1)
        double qsecond = 2*b*b - 8*a*b + 2*a*a;
        return s.evalSecondDer(x) + e[i]*qsecond;
    }
};

static void writeSVG(const string& filename,
                      const vector<pair<double,double>>& splinePts,
                      const vector<pair<double,double>>& splineDerPts,
                      const vector<Point>& dataPts,
                      const vector<pair<double,double>>* integralPts,
                      double xmin, double xmax, double ymin, double ymax)
{
    int W=900, H=600;

    auto X = [&](double x){ return (x - xmin)/(xmax-xmin)*W; };
    auto Y = [&](double y){ return H - (y - ymin)/(ymax-ymin)*H; };

    ofstream out(filename);
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << W << "\" height=\"" << H << "\">\n";
    out << "<rect x=\"0\" y=\"0\" width=\"" << W << "\" height=\"" << H << "\" fill=\"white\"/>\n";

    // Axes
    out << "<line x1=\"0\" y1=\"" << Y(0) << "\" x2=\"" << W << "\" y2=\"" << Y(0)
        << "\" stroke=\"#ddd\"/>\n";

    // Spline curve
    auto toPoints = [&](const vector<pair<double,double>>& poly){
        ostringstream oss;
        for(size_t i=0;i<poly.size();i++){
            oss << (int)llround(X(poly[i].first)) << "," << (int)llround(Y(poly[i].second));
            if(i+1<poly.size()) oss << " ";
        }
        return oss.str();
    };

    out << "<polyline fill=\"none\" stroke=\"#1f77b4\" stroke-width=\"2\" points=\""
        << toPoints(splinePts) << "\"/>\n";

    // Derivative curve (optional visualization)
    out << "<polyline fill=\"none\" stroke=\"#d62728\" stroke-width=\"1.5\" stroke-dasharray=\"5,4\" points=\""
        << toPoints(splineDerPts) << "\"/>\n";

    // Integral curve (optional)
    if(integralPts){
        out << "<polyline fill=\"none\" stroke=\"#2ca02c\" stroke-width=\"2\" points=\""
            << toPoints(*integralPts) << "\"/>\n";
    }

    // Data points
    for(const auto& p : dataPts){
        out << "<circle cx=\"" << (int)llround(X(p.x)) << "\" cy=\"" << (int)llround(Y(p.y))
            << "\" r=\"4.5\" fill=\"#000\" stroke=\"#fff\" stroke-width=\"1.2\"/>\n";
    }

    // Tangent lines at data points: y = yp_i*(x-x_i)+y_i
    // draw short segment in x around each knot
    double span = xmax - xmin;
    double dxSeg = 0.03 * span; // segment length in x-direction
    for(const auto& p : dataPts){
        double xA = max(xmin, p.x - dxSeg);
        double xB = min(xmax, p.x + dxSeg);

        double yA = p.y + p.yp*(xA - p.x);
        double yB = p.y + p.yp*(xB - p.x);

        out << "<line x1=\"" << (int)llround(X(xA)) << "\" y1=\"" << (int)llround(Y(yA))
            << "\" x2=\"" << (int)llround(X(xB)) << "\" y2=\"" << (int)llround(Y(yB))
            << "\" stroke=\"#ff7f0e\" stroke-width=\"1.3\"/>\n";
    }

    out << "</svg>\n";
}

int main(int argc, char** argv){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string inputFile = "data.txt";
    string outFile = "out.svg";
    string mode = "cubic"; // cubic | quartic
    bool drawIntegral = false;

    for(int i=1;i<argc;i++){
        string a=argv[i];
        if(a=="--in" && i+1<argc) inputFile=argv[++i];
        else if(a=="--out" && i+1<argc) outFile=argv[++i];
        else if(a=="--mode" && i+1<argc) mode=argv[++i];
        else if(a=="--integral") drawIntegral=true;
    }

    ifstream in(inputFile);
    if(!in){
        cerr << "Failed to open " << inputFile << "\n";
        return 1;
    }

    int n; in >> n;
    vector<Point> pts(n);
    for(int i=0;i<n;i++){
        in >> pts[i].x >> pts[i].y >> pts[i].yp;
    }

    CubicSpline cubic;
    QuarticSpline quartic;

    if(mode=="quartic"){
        quartic.build(pts);
    } else {
        cubic.build(pts);
    }

    double xmin0=pts.front().x, xmax0=pts.back().x;

    // sample points
    int N=1200;
    vector<pair<double,double>> splineY, splineYp;
    splineY.reserve(N);
    splineYp.reserve(N);

    vector<pair<double,double>> integralPts;

    vector<double> xs(N), ys(N), yps(N);
    for(int k=0;k<N;k++){
        double x = xmin0 + (xmax0-xmin0)*k/(N-1);
        double y, yp;
        if(mode=="quartic"){
            y = quartic.eval(x);
            yp = quartic.evalDer(x);
        } else {
            y = cubic.eval(x);
            yp = cubic.evalDer(x);
        }
        xs[k]=x; ys[k]=y; yps[k]=yp;
    }

    if(drawIntegral){
        vector<double> integral(N,0.0);
        for(int k=1;k<N;k++){
            double dx = xs[k]-xs[k-1];
            integral[k] = integral[k-1] + 0.5*(ys[k]+ys[k-1])*dx;
        }
        integralPts.reserve(N);
        for(int k=0;k<N;k++) integralPts.push_back({xs[k], integral[k]});
    }

    for(int k=0;k<N;k++){
        splineY.push_back({xs[k], ys[k]});
        splineYp.push_back({xs[k], yps[k]});
    }

    // Compute bounds including tangent lines and integral (if shown)
    double xmin=xmin0, xmax=xmax0, ymin=1e300, ymax=-1e300;

    // Include spline
    for(auto &pr: splineY){ ymin=min(ymin, pr.second); ymax=max(ymax, pr.second); }
    // Include data points
    for(auto &p: pts){ ymin=min(ymin, p.y); ymax=max(ymax, p.y); }

    // Include tangents (approx: evaluate at a couple x offsets)
    double span = xmax0 - xmin0;
    double dxSeg = 0.03 * span;
    for(auto &p: pts){
        double xA = max(xmin0, p.x - dxSeg);
        double xB = min(xmax0, p.x + dxSeg);
        double yA = p.y + p.yp*(xA - p.x);
        double yB = p.y + p.yp*(xB - p.x);
        ymin=min(ymin, min(yA,yB));
        ymax=max(ymax, max(yA,yB));
    }

    // If derivative curve included in same y-scale, expand bounds a bit for it
    for(auto &pr: splineYp){ ymin=min(ymin, pr.second); ymax=max(ymax, pr.second); }

    // If integral included, also expand bounds
    if(drawIntegral){
        for(auto &pr: integralPts){ ymin=min(ymin, pr.second); ymax=max(ymax, pr.second); }
    }

    if(!(ymax>ymin)) { ymax=ymin+1; }
    double padX = 0.05*(xmax-xmin);
    double padY = 0.10*(ymax-ymin);
    xmin -= padX; xmax += padX; ymin -= padY; ymax += padY;

    writeSVG(outFile, splineY, splineYp, pts, drawIntegral ? &integralPts : nullptr,
             xmin, xmax, ymin, ymax);

    cerr << "Wrote " << outFile << "\n";
    return 0;
}