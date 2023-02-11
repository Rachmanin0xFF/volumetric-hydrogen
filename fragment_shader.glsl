// @author Adam Lastowka

//================= MATH =================//

precision highp float;

#define PI 3.14159265359
#define TWO_PI 6.28318530718
// The maximum order of any polynomial used to calculate the wavefunction
#define ORDER 24

// I represent polynomials as a list of coefficients.
// e.g. y(x) = c[0] + c[1]*x + c[2]*x^2 + ... + c[n]*x^n
struct Polynomial
{
    float c[ORDER];
};

Polynomial derivative(Polynomial p) {
    Polynomial o;
    for(int i = 1; i < ORDER; i++) {
        o.c[i-1] = p.c[i]*float(i);
    }
    return o;
}
Polynomial multx(Polynomial p) {
    Polynomial o;
    for(int i = 0; i < ORDER-1; i++) {
        o.c[i+1] = p.c[i];
    }
    return o;
}
Polynomial mult(Polynomial p, float f) {
    Polynomial o;
    for(int i = 0; i < ORDER; i++) {
        o.c[i] = p.c[i]*f;
    }
    return o;
}
Polynomial add(Polynomial a, Polynomial b) {
    Polynomial o;
    for(int i = 0; i < ORDER; i++) {
        o.c[i] = a.c[i] + b.c[i];
    }
    return o;
}
float eval(Polynomial p, float x) {
    float sum = p.c[0];
    float xx = x;
    for(int i = 1; i < ORDER; i++) {
        sum += p.c[i]*xx;
        xx *= x;
    }
    return sum;
}
// returns (-1)^n
int m1pow(int n) {
    return n%2==0?1:-1;
}

struct Orbital {
    Polynomial P;
    Polynomial L;
    int n;
    int l;
    int m;
};

// Many thanks to Justin Willmert for his writeup:
// https://justinwillmert.com/articles/2020/pre-normalizing-legendre-polynomials/

// Note: this function does NOT contain the (1-x^2)^(m/2) term; this factor is added
// later in the eval_orbital() function

Polynomial associated_Legendre(int l, int m) {
    Polynomial o;
    if(m > l) return o;
    o.c[0] = 0.28209479177; // SH NORM = sqrt(1/(4PI))
    
    int ll = 0;
    int mm = 0;
    for(; mm < m; mm++) {
        ll++;
        float u_l = sqrt(1.0+1.0/float(2*ll));
        o = mult(o, -u_l);
    }
    if(l==m) return o;
    
    Polynomial m1 = o;
    ll++;
    float v_l = sqrt(float(1 + 2*ll)); // SH NORM
    o = mult(multx(o), v_l);
    if(l==ll) return o;
    
    Polynomial m0 = o;
    for(;ll < l;) {
        ll++;
        float c = float(2*ll+1)/(float((2*ll-3)*(ll*ll-mm*mm)));
        float a = sqrt(c*float(4*(ll-1)*(ll-1)-1)); // SH NORM
        float b = sqrt(c*float((ll-1)*(ll-1)-mm*mm)); // SH NORM
        
        o = add(mult(multx(m0), a), mult(m1, -b));
        m1 = m0;
        m0 = o;
    }
    return o;
}

// I tried Willmert's prenormalization technique for these, too,
// but the expressions for the normalizing factors were so hefty that I don't think
// they would provide any significant advantage (computationally) over this method.
//
// I just used a CAS -- it might still be worth exploring later.

Polynomial generalized_Laguerre(int n, int l) {
    Polynomial o;
    o.c[0] = 1.0;
    if(n-l-1<=0) return o;
    
    Polynomial m1 = o;
    float a = float(2*l+1);
    o.c[0] = float(1.0+a); o.c[1] = -1.0;
    
    Polynomial m0 = o;
    for(int kk = 1; kk<(n-l-1); kk++) {
        o = add(mult(multx(m0), -1.0), mult(m0, float(2*kk+1)+a));
        o = add(o, mult(m1, -float(kk)-a));
        o = mult(o, 1.0/float(kk+1));
        m1 = m0;
        m0 = o;
    }
    return o;
}

// Is this accurate?
// No, definitely not; (psi*)(psi) does not integrate to 1 with this factor.
// But the accurate normalization factor leads to some orbitals looking much dimmer.
// Rather than correcting this later, I just use a different "wrong" normalization.

float radial_normalization(Orbital o) {
    float N = float(o.n*o.n)*0.125;
    for(int i = o.n-o.l; i <= o.n+o.l; i++) N *= float(i);
    return sqrt(1.0/N);
}

Orbital get_orbital(int n, int l, int m) {
    Orbital o;
    o.n = n;
    o.l = l;
    o.m = m;
    o.P = associated_Legendre(l, m);
    o.L = mult(generalized_Laguerre(n, l), radial_normalization(o));
    return o;
}

float eval_orbital_radial(Orbital o, float r) {
    float rr = 2.0*r/float(o.n);
    float L = eval(o.L, rr);
    L *= pow(rr, float(o.l));
    L *= exp(-r/float(o.n));
    return L;
}

// sph_coords = (theta, phi)
// theta is the polar angle
// phi is the azimuthal angle

vec2 eval_orbital_SH(Orbital o, vec2 sph_coords) {
    float x = cos(sph_coords.x);
    float r = eval(o.P, x)*pow(1.0-x*x, float(o.m)*0.5);
    return r*vec2(cos(float(o.m)*sph_coords.y), sin(float(o.m)*sph_coords.y));
}

vec2 eval_orbital(Orbital o, vec3 pos) {
    float R = length(pos.xy);
    float rho = length(pos);
    vec2 sph_coords = vec2(atan(R, pos.z), atan(pos.y, pos.x));
    
    vec2 Y = eval_orbital_SH(o, sph_coords);
    float radial = eval_orbital_radial(o, rho);
    
    return Y*radial;
}


//================= COLOR =================//

// You can find these and more here:
// https://github.com/Rachmanin0xFF/GLSL-Color-Functions
// (MIT license, 2022)

//                          0.3127/0.3290  1.0  (1.0-0.3127-0.3290)/0.329
const vec3 D65_WHITE = vec3(0.95045592705, 1.0, 1.08905775076);

vec3 WHITE = D65_WHITE;

// Chromatic adaptation between D65<->D50
// XYZ color space does not depend on a reference white, but all other matrices here
// assume D65. These "restretch" XYZ to the D50 reference white so the others can sitll work with D50.

// from https://www.color.org/sRGB.pdf
const mat3 XYZ_TO_XYZ50_M = mat3(
    1.0479298208405488, 0.022946793341019088, -0.05019222954313557,
    0.029627815688159344, 0.990434484573249, -0.01707382502938514,
    -0.009243058152591178, 0.015055144896577895, 0.7518742899580008
);
const mat3 XYZ50_TO_XYZ_M = mat3(
    0.9554734527042182, -0.023098536874261423, 0.0632593086610217,
    -0.028369706963208136, 1.0099954580058226, 0.021041398966943008,
    0.012314001688319899, -0.020507696433477912, 1.3303659366080753
);

// RGB<->XYZ
// from IEC 61966-2-1:1999/AMD1:2003 (sRGB color amendment 1)
const mat3 RGB_TO_XYZ_M = mat3(
    0.4124, 0.3576, 0.1805,
    0.2126, 0.7152, 0.0722,
    0.0193, 0.1192, 0.9505
);
const mat3 XYZ_TO_RGB_M = mat3(
    3.2406255, -1.5372080, -0.4986286,
    -0.9689307, 1.8757561, 0.0415175,
    0.0557101, -0.2040211, 1.0569959
);

// sRGB<->RGB
// sRGB is standard "monitor" space, and the standard colorspace of the internet.
// The EOTF is roughly equivalent to a gamma of 2.2, but it acts differently in low values.
float UNCOMPAND_SRGB(float a) {
    return (a > 0.04045) ? pow((a + 0.055) / 1.055, 2.4) : (a / 12.92);
}
vec3 SRGB_TO_RGB(vec3 srgb) {
    return vec3(UNCOMPAND_SRGB(srgb.x), UNCOMPAND_SRGB(srgb.y), UNCOMPAND_SRGB(srgb.z));
}
float COMPAND_RGB(float a) {
    return (a <= 0.0031308) ? (12.92 * a) : (1.055 * pow(a, 0.41666666666) - 0.055);
}
vec3 RGB_TO_SRGB(vec3 rgb) {
    return vec3(COMPAND_RGB(rgb.x), COMPAND_RGB(rgb.y), COMPAND_RGB(rgb.z));
}

// RGB<->XYZ
// XYZ is the classic tristimulus color space developed in 1931 by the International Commission on Illumination (CIE, confusingly).
// Most conversions between color spaces end up going through XYZ; it is a central 'hub' in the color space landscape.
vec3 RGB_TO_XYZ(vec3 rgb) {
    return WHITE == D65_WHITE ? (rgb * RGB_TO_XYZ_M) : ((rgb * RGB_TO_XYZ_M) * XYZ_TO_XYZ50_M);
}
vec3 XYZ_TO_RGB(vec3 xyz) {
    return WHITE == D65_WHITE ? (xyz * XYZ_TO_RGB_M) : ((xyz * XYZ50_TO_XYZ_M) * XYZ_TO_RGB_M);
}

// L*a*b*/CIELAB
// CIELAB was developed in 1976 in an attempt to make a perceptually uniform color space.
// While it doesn't always do a great job of this (especially in the deep blues), it is still frequently used.
float XYZ_TO_LAB_F(float x) {
    //          (24/116)^3                         1/(3*(6/29)^2)     4/29
    return x > 0.00885645167 ? pow(x, 0.333333333) : 7.78703703704 * x + 0.13793103448;
}
vec3 XYZ_TO_LAB(vec3 xyz) {
    vec3 xyz_scaled = xyz / WHITE;
    xyz_scaled = vec3(
        XYZ_TO_LAB_F(xyz_scaled.x),
        XYZ_TO_LAB_F(xyz_scaled.y),
        XYZ_TO_LAB_F(xyz_scaled.z)
    );
    return vec3(
        (116.0 * xyz_scaled.y) - 16.0,
        500.0 * (xyz_scaled.x - xyz_scaled.y),
        200.0 * (xyz_scaled.y - xyz_scaled.z)
    );
}
float LAB_TO_XYZ_F(float x) {
    //                                     3*(6/29)^2         4/29
    return (x > 0.206897) ? x * x * x : (0.12841854934 * (x - 0.137931034));
}
vec3 LAB_TO_XYZ(vec3 Lab) {
    float w = (Lab.x + 16.0) / 116.0;
    return WHITE * vec3(
        LAB_TO_XYZ_F(w + Lab.y / 500.0),
        LAB_TO_XYZ_F(w),
        LAB_TO_XYZ_F(w - Lab.z / 200.0)
    );
}

// LCh
// LCh is simply L*a*b* converted to polar coordinates.
// Note: by convention, h is in degrees!
vec3 LAB_TO_LCH(vec3 Lab) {
    return vec3(
        Lab.x,
        sqrt(dot(Lab.yz, Lab.yz)),
        atan(Lab.z, Lab.y) * 57.2957795131
    );
}
vec3 LCH_TO_LAB(vec3 LCh) {
    return vec3(
        LCh.x,
        LCh.y * cos(LCh.z * 0.01745329251),
        LCh.y * sin(LCh.z * 0.01745329251)
    );
}

// Composite function one-liners
vec3 SRGB_TO_XYZ(vec3 srgb) { return RGB_TO_XYZ(SRGB_TO_RGB(srgb)); }
vec3 XYZ_TO_SRGB(vec3 xyz)  { return RGB_TO_SRGB(XYZ_TO_RGB(xyz));  }

vec3 SRGB_TO_LAB(vec3 srgb) { return XYZ_TO_LAB(SRGB_TO_XYZ(srgb)); }
vec3 LAB_TO_SRGB(vec3 lab)  { return XYZ_TO_SRGB(LAB_TO_XYZ(lab));  }

vec3 SRGB_TO_LCH(vec3 srgb) { return LAB_TO_LCH(SRGB_TO_LAB(srgb)); }
vec3 LCH_TO_SRGB(vec3 lch)  { return LAB_TO_SRGB(LCH_TO_LAB(lch));  }


//================= RAY MARCHING =================//


// RGB for color, alpha for density
vec4 worldVolumetric(vec3 pos, Orbital o) {
    pos *= 13.0;
    float theta = atan(pos.x, pos.z);
    float phi = atan(length(pos.xz), pos.y);
    vec2 sh = vec2(theta, phi)/(length(pos)+1.0)/10.0;
    
    sh = eval_orbital(o, pos);
    
    vec3 col = XYZ_TO_RGB(LAB_TO_XYZ(LCH_TO_LAB(vec3(75.0, 52.0, atan(sh.y, sh.x)/PI*180.0+90.0))));
 
    //sh.y *= 0.0;
    float lcol = 1.0*dot(sh, sh);
    //if(pos.x > 0.0) return vec4(0.0);
    return vec4(col*0.7, 3000.0*lcol);
}

struct RaymarchResult {
    vec3 position;
    vec3 color;
    float traversed;
};

float march_shadow(in vec3 ro, in vec3 rd, int max_iter) {
    float step_size = 0.1;
    float amt = 0.0;
    ro += rd*(1.0 + texture(iChannel1, rd.xy).r);
    int i = 0;
    for(i = 0; i < max_iter; i++) {
        vec4 samp = worldVolumetric(ro, get_orbital(3, 3, 3));
        amt += samp.a*step_size;
    	ro += rd*step_size;
    }
    return exp(-amt);
}
    
RaymarchResult march(in vec3 ro, in vec3 rd, int max_iter, in vec2 fc) {
    Orbital orb = get_orbital(5, 3, 2);
    vec4 col = vec4(0.0);
    RaymarchResult data;
    
    float step_size = 0.2;
    
    vec3 light_dir = vec3(-1.0, 0.1, -0.3);
    
    ro += rd*(1.0 + step_size*texture(iChannel1, fc.xy/iChannelResolution[1].xy).r);
    int i = 0;
    for(i = 0; i < max_iter; i++) {
        vec4 samp = worldVolumetric(ro, orb);
        samp.rgb = max(samp.rgb, vec3(0.0));
        
        float shadow = 1.0;
        data.traversed += samp.a*step_size;
        
        col += samp*step_size*samp.a*(1.0 - col.a);
        
        if(col.a > 1.0) break;
    	ro += rd*step_size;
    }
    data.position = ro;
    data.color = col.rgb;
    return data;
}

//================= CAMERA =================//

struct Camera
{
	vec3 position;

	// Camera space stuff
	vec3 forwards;
	vec3 left; // remove left/up to save space in future
	vec3 up;

	vec3 rayDir;
};

Camera getCam()
{
	Camera cam;
	vec3 lookAt = vec3(0.0, 0.0, 0.0);
	cam.position = vec3(7.0*cos(iTime), sin(iTime)*3.0, 7.0*sin(iTime));

	// figure out camera space from position and lookAt
	cam.up = vec3(0, 1, 0);
	cam.forwards = normalize(lookAt - cam.position);
	cam.left = normalize(cross(cam.forwards, cam.up));
	cam.up = normalize(cross(cam.left, cam.forwards));

	// find view ray - fustrum intersection for this pixel
	vec3 fustrumFront = cam.position + cam.forwards;
	vec2 screenSpace = (2.0*gl_FragCoord.xy/iResolution.xy - 1.0)*0.45;
	float aspect = iResolution.x/iResolution.y;
	vec3 fustrumIntersect = fustrumFront + screenSpace.x*cam.left*aspect + screenSpace.y*cam.up;

	// direction to march in
	cam.rayDir = normalize(fustrumIntersect-cam.position);

	return cam;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    Camera cam = getCam();
    RaymarchResult rm = march(cam.position, cam.rayDir, 60, fragCoord);
    fragColor.xyz = RGB_TO_SRGB(rm.color);
}
