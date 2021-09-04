		uniform vec3 iResolution;
		uniform float iTime;
        uniform sampler2D buffer;
        uniform sampler2D texture01;
        uniform sampler2D texture02;
        uniform sampler2D texture03;
        uniform sampler2D texture04;
        uniform sampler2D texture05;
        uniform sampler2D texture06;

        uniform vec3 cameraPos;
        uniform vec3 cameraDir;
        uniform vec3 cameraLens;
        uniform bool LensCheck;

        uniform vec3 _BASECOLOR;
        uniform float _SUBSURFACE;
        uniform float _METALLIC;
        uniform float _SPECULAR;
        uniform float _SPECULARTINT;
        uniform float _ROUGHNESS;
        uniform float _ANISTROPIC;
        uniform float _SHEEN;
        uniform float _SHEENTINT;
        uniform float _CLEARCOAT;
        uniform float _CLEARCOATTINT;
        uniform vec3 _POSITION;
        uniform float _LUMINANCE;
        uniform float _SPHERE;

        uniform int World;
        uniform float WorldLumi;
        uniform bool Normaltest;
#define M_PI float(3.141592)
#define HASHSCALE3 vec3(.1031, .1030, .0973)
#define HASHSCALE2 vec2(.1392,.1953)

//Rayや球、衝突情報の構造体の宣言
struct Ray{
    vec3 origin;
    vec3 direction;
};
struct Sphere{
    float radius;
    vec3 position;
    vec3 color;
    int number;
};
struct Info{
    float distance;
    vec3 normal;
    vec3 position;
    vec3 color;
    vec2 uv;
    int number;
};

//線形補完
vec3 lerp(vec3 v1,vec3 v2,float t){
    return v1 + t * (v2 - v1);
}
float lerp(float a, float b, float t){
    return a + t * (b - a);
}
//ランダム関数、サンプリングなどに使う。
float random (vec2 st) {
    return fract(dot(st.xy,
                         vec2(12.9898,78.233)));
}

//reffer https://github.com/gam0022/webgl/blob/master/pathtracing_sandbox.html by gam0022
vec2 hash22(vec2 p2){
    p2 = fract(p2 * HASHSCALE2);
    p2 += dot(p2,p2.yx + 12.42);
    return fract((p2.xx + p2.yx) * p2);
}
vec2 hash23(vec3 p3)
{
	p3 = fract(p3 * HASHSCALE3);
	p3 += dot(p3, p3.yzx+19.19);
	return fract((p3.xx+p3.yz)*p3.zy);
}
//接空間の基底を求める関数
void tangentSpaceBasis(vec3 normal,inout vec3 t,inout vec3 b){
    if (abs(normal.y) < 0.9)
    {
        t = cross(normal, vec3(0, 1, 0));
    }
    else
    {
        t = cross(normal, vec3(0, 0, -1));
    }
    t = normalize(t);
    b = cross(t, normal);
    b = normalize(b);
}

//球のＵＶを取得する
vec2 sphereUV(vec3 direction){
    vec3 pos = normalize(direction);
    float theta = acos(pos.y);
    float phi = acos(pos.x / sqrt(pos.z * pos.z + pos.x * pos.x));
    if(pos.z < 0.0){
    phi *= -1.;
    }
    float u = clamp(theta / M_PI,0.,1.);
    float v = (phi + M_PI) / (2.0 * M_PI);
    v = clamp(1.0 - v,0.,1.0);
    u = clamp(1.0 - u,0.,1.0);
    return vec2(v,u);
}

//ワールド座標からローカル座標への変換
vec3 worldtoLoacal(vec3 v,vec3 lx, vec3 ly,vec3 lz){
    return vec3(v.x * lx.x + v.y* lx.y + v.z * lx.z,
                 v.x * ly.x + v.y * ly.y + v.z * ly.z,
                 v.x * lz.x + v.y * lz.y + v.z * lz.z);
}

//ローカルからワールド座標への変換
vec3 localToWorld(const vec3 v, const vec3 lx, const vec3 ly,
                   const vec3 lz)
{
    return vec3(v.x * lx.x + v.y * ly.x + v.z * lz.x,
                 v.x * lx.y + v.y * ly.y + v.z * lz.y,
                 v.x * lx.z + v.y * ly.z + v.z * lz.z);
}

//コサイン重点サンプリング
vec3 sampleHemisphere(float u, float v, inout float pdf){
	float theta = acos(clamp(1.0 - 2.0 * u, -1.0, 1.0)) / 2.0;
    	float phi = 2.0 * M_PI * v;
    	pdf = cos(theta) / M_PI;
    	return vec3(cos(phi) * sin(theta),cos(theta),
                 sin(phi) * sin(theta));
}

//https://hal.archives-ouvertes.fr/hal-01509746/document
vec3 sampleGGXheitz(float u,float v, vec3 wo,inout vec3 normal){
    vec3 T1 = (wo.z <0.99) ? normalize(cross(wo,vec3(0.0,0.0,1.0))) : vec3(1.0,0.0,0.0);
    vec3 T2 = cross(T1,wo);

    float a = 1.0 / (1.0 + wo.z);
    float r = sqrt(u);
    float phi = (v < a) ? v / a * M_PI : M_PI + (v - a) / (1.0 - a) * M_PI;
    float P1 = r * cos(phi);
    float P2 = r * sin(phi) * ((v < a) ? 1.0 : wo.z);
    
    vec3 P = P1 * T1 + P2 * T2;
    vec3 N = P + sqrt(max(0.0, 1.0 - length(P) * length(P))) * wo;
    N = normalize(N);
    normal =N;
    vec3 wi = reflect(wo,N);
    return wi;
}

vec3 sampleGGXwalter(float u,float v ,float alpha,vec3 wo){
    float theta = atan(alpha * sqrt(u) / sqrt(max(0.0,1.0 - u)));
    float phi = 2. * M_PI * v;
    vec3 normal = vec3(cos(phi) * sin(theta),cos(theta),
                 sin(phi) * sin(theta));
    return normalize(reflect(wo,normal));
}

//球の衝突判定
bool intersectSphere(Ray R,Sphere S,inout Info info){
    vec3 a = R.origin - S.position;
    float b = dot(a,R.direction);
    float c = dot(a,a) - (S.radius * S.radius);
    float d = b * b - c;
    if( d <= 0.0){
        return false;
    }

    float t1 = -b - sqrt(d);
    float t2 = -b + sqrt(d);

    float t = t1;
    if(t < 0.00001 || t > 10000.){
        t = t2;
        if(t < 0.0001 || t > 10000.){
            return false;
        } 
    }

    //if(t > info.distance) return false;

    info.distance = t;
    info.position = R.origin + t * R.direction;
    info.normal = normalize(info.position - S.position);
    info.color = vec3(1.0);
    info.number = S.number;
    info.uv = sphereUV(info.normal);
    return true;
}

//平面の衝突判定
bool cercleplaneIntersect(Ray ray,inout Info info,float radiance){
    vec3 normal = normalize(vec3(0.,1.,0.));
    
    float sn = dot(ray.origin,normal);
    float dn = dot(ray.direction,normal);

    if(abs(dn) < 0.001) return false;

    float t = - sn/ dn;
    if(t < 0.0001 || t > 10000.) return false;
    vec3 pos = ray.origin + t * ray.direction;
    if(length(pos) > radiance) return false;
    info.distance = t;
    info.normal = normal;
    if(dot(normal,ray.direction) > 0.0) info.normal = -normal;
    info.position = pos;
    info.color = vec3(0.9);
    info.number = 0;
    return true;
}
//四角い平面の衝突判定
bool squarePlaneIntersect(Ray ray, inout Info info){
    vec3 normal = normalize(vec3(0.,0.,1.));
    vec3 PlanePos = vec3 (0,3,-3);
    float rn = dot(normal,PlanePos);
    float sn = dot(ray.origin,normal);
    float dn = dot(ray.direction,normal);

    if(abs(dn) < 0.001) return false;
    float t = (rn - sn) / dn;
    if(t < 0.0001 || t > 10000.) return false;

    vec3 pos = ray.origin + t * ray.direction;
    vec3 lim = pos - PlanePos;
    float limx = 3.5,limy = 2.5;
    if(abs(lim.y) > limy || abs(lim.x) > limx) return false;

    info.distance = t;
    info.normal = normal;
    if(dot(normal,ray.direction) > 0.0) info.normal = -normal;
    info.color = vec3(0.0);
    info.position = pos;
    info.number = 10;
    info.uv = vec2( 0.5 * (limx + lim.x)/limx,0.5 *(limy + lim.y) / limy);
    return true;
}
vec3 WorldSphere(Ray ray){
    vec3 a = ray.origin;
    float b = dot(a,ray.direction);
    float c = dot(a,a) - (10000. * 10000.);
    float d = b * b - c;
    if( d <= 0.0){
        return vec3(0.);
    }

    float t1 = -b - sqrt(d);
    float t2 = -b + sqrt(d);

    float t = t1;
    if(t < 0.0001 || t > 100000.){
        t = t2;
        if(t < 0.0001 || t > 100000.){
            return vec3(0.);
        } 
    }
    vec3 hitpos = ray.origin + ray.direction * t;
    vec2 worlduv = sphereUV(hitpos);
    vec3 radiance = texture2D(texture05,worlduv).xyz * WorldLumi;
    if(World == 1){
        radiance = vec3(0.0);
    }
    else if(World == 2){
        radiance = vec3(0.9);
    }else if(World == 3){
        radiance = texture2D(texture06,worlduv).xyz;
    }
    return radiance;
}

//シーンの衝突判定
bool Scene(Ray r,inout Info info){

    Sphere sphere;
    sphere.radius = 1.0;
    sphere.position = vec3(0.0,1.01,0.0);
    sphere.number = 1;
    float width = 6.0;
    float spherenumber = 5.;
    float dis = width * 2.0 / (spherenumber-1.);    
    bool hit = false;
    info.distance = 10000.;
    for(int i = 1; i <= 5; i++){
        sphere.position = vec3(-width + dis * float(i-1),1.0,0.0);
        sphere.number = i;
        Info tes;
        tes.distance = 10000.;
        if(intersectSphere(r,sphere,tes)){
            if(tes.distance < info.distance){
                hit = true;
                info = tes;
            }
        }
    }

    Info tes;
    tes.distance = 10000.;
    if(cercleplaneIntersect(r,tes,8.0)){
        if(tes.distance < info.distance){
            hit = true;
            info = tes;
        }
    }

    Info test;
    test.distance = 10000.;
    sphere.radius = _SPHERE;
    sphere.number = 6;
    sphere.position = _POSITION;
    if(intersectSphere(r,sphere,test)){
        if(test.distance <info.distance){
            hit = true;
            info = test;
        }
    }

    //Info test;
    //test.distance = 10000.;
    //if(squarePlaneIntersect(r,test)){
    //    if(test.distance < info.distance){
    //        hit = true;
    //        info = test;
    //    }
    //}
    
    //Info l;
    //sphere.position = _POSITION;
    //sphere.number = 6;
    //l.distance = 10000.;
    //if(intersectSphere(r,sphere,l)){
    //    if(l.distance < info.distance){
    //        hit = true;
    //        info = l;
    //    }
    //}
    return hit;
}

//BRDF
//Lamvert
vec3 Lambert(vec3 rho,inout vec3 wi,inout float pdf,vec2 random){
     wi = sampleHemisphere(random.x,random.y,pdf);
     return rho / M_PI;
}
// Perfect Reflect
vec3 reflect(vec3 wo,vec3 n){
    return normalize(-wo + 2. * dot(wo,n) * n);
}

vec3 Mirror(vec3 rho,vec3 wo,inout vec3 wi,inout float pdf){
    wi = reflect(wo,vec3(0.,1.0,0.));
    pdf = 1.0;
    return rho / abs(wi[1]);
}

// Non Roughness Glass
float length2(vec3 v){
    return length(v) * length(v);
}
vec3 refract(vec3 wo,vec3 n,float n1, float n2){
    vec3 th = - (n1 / n2) * (wo - dot(wo,n) * n);
    if(length2(th)> 1.0){
        return reflect(wo,n);
    }
    vec3 tp = -sqrt(max(1.0f - length2(th),0.0)) * n;
    return normalize(th + tp);
}
float glassFresnel(vec3 wo,vec3 n,float ior1,float ior2){
    float f0 = pow((ior1 - ior2) / (ior1 + ior2),2.0);
    return f0 + (1.0f - f0) * pow(1.0 - abs(dot(wo,n)),5.0);
}

vec3 SmoothGlass(vec3 wo,inout vec3 wi,inout float pdf,float random){
    vec3 rho = vec3(0.9);
    float ior1,ior2;
    float glassior = 1.33;
    vec3 n = vec3(0.,1.,0.);
    if(wo.y > 0.0){
        ior1 = 1.0;
        ior2 = glassior = 1.33; 
        n = vec3(0.,1.,0.);
    }
    else
    {
        ior1 = glassior;
        ior2 = 1.0;
        n = vec3(0.,-1.,0.);
    }

    float fr = glassFresnel(wo,n,ior1,ior2);
    if(random < fr){
        wi = reflect(wo,n);
        pdf = 1.0;
        return rho / abs(wi.y);
    }
    else{
        pdf = 1.0;
        wi = refract(wo,n,ior1,ior2);
        return rho / abs(wi.y);
    }
}


//MicrofasetBRDF
// https://qiita.com/aa_debdeb/items/f813bdcbd8524a66a11b
vec3 fresnelSchlick(vec3 f0,float cosine){
    return f0 + (1.0 - f0) * pow(1.0 - cosine,5.0);
}

float DGGX(float nh,float alpha){
    float d = (1.0 - (1.0 - alpha * alpha) * nh * nh);
    d = max(d,0.01);
    return alpha * alpha /(M_PI * d * d);
}

float Lambda(float alpha,float xn){
    return 0.5 * (-1.0 + sqrt(1.0 + alpha * alpha *(1./(xn * xn) - 1.0)));
}

float GSJMSF(float nv,float nl,float alpha){
    float lambdav = Lambda(alpha,nv);
    float lambdal = Lambda(alpha,nl);
    return 1./(1. + lambdav + lambdal);
}

vec3 MicroFacet(vec3 F0, vec3 wo,inout vec3 wi,float roughness,inout float pdf,vec2 random){
    vec3 n = vec3(0.0,1.0,0.0);
    vec3 N;
    wi = sampleHemisphere(random.x,random.y,pdf);
    //wi = sampleGGXheitz(random.x,random.y,wo,N);    
    float alpha = roughness * roughness;
    //wi = sampleGGXwalter(random.x,random.y,alpha,wo);
    vec3 v = wo;
    vec3 l = wi;
    vec3 h = normalize(v+l);
    float nh = dot(n,h);
    float nl = dot(n,l);
    float nv = dot(n,v);
    float vh = dot(v,h);
    
    float D = DGGX(nh,alpha);
    float G = GSJMSF(nv,nl,alpha);
    vec3 F = fresnelSchlick(F0,vh);
    
    vec3 f = F* D * G / (4.0 * nv * nl );
    //pdf = D * nh / (4.0 * vh);
    //pdf = nh / (4.0 * vh);
    return f;
}

//Disney BRDF
float SchlickFresnel(float w,float F90){
    float m = clamp(max(1.-w,0.0),0.,1.);
    float wn5 = m * m * m * m * m;
    return 1.0 + (F90 - 1.0) * wn5;
}

vec3 Spec(vec3 wo,vec3 wi,float roughness,vec3 F0){
    vec3 n = vec3(0.0,1.0,0.0);
    vec3 v = wo;
    vec3 l = wi;
    vec3 h = normalize((v+l)/2.0);
    float nh = dot(n,h);
    float nl = dot(n,l);
    float nv = dot(n,v);
    float vh = dot(v,h);
    
    float alpha = roughness;
    
    float D = DGGX(nh,alpha);
    float G = GSJMSF(nv,nl,alpha);
    vec3 F = fresnelSchlick(F0,vh);
    
    vec3 f = F* D * G / (4.0 * nv * nl );
    return f;
}
//Berry
float DBerry(float NdotH,float a){
    float a2 = a*a;
    float t = 1.0 + (a2 - 1.0) * NdotH * NdotH;
    return (a2 - 1.0) / (M_PI * log(a2) * t);
}

vec3 clearC(vec3 wo,vec3 wi,float alpha,float cl){
    vec3 n = vec3(0.0,1.0,0.0);
    vec3 v = wo;
    vec3 l = wi;
    vec3 h = normalize((v+l)/2.0);
    float nh = dot(n,h);
    float nl = dot(n,l);
    float nv = dot(n,v);
    float vh = dot(v,h);
    
    float D = DBerry(nh,alpha);
    float G = GSJMSF(nv,nl,0.25);
    vec3 F = fresnelSchlick(vec3(0.04),vh);
    
    vec3 f = 0.25 * cl * F* D * G / (4.0 * nv * nl );
    return f;
}
vec3 DisneyBRDF(vec3 wo,inout vec3 wi,vec2 random,inout float pdf,vec3 baseColor,float subsurface,float metallic,float specular,float specularTint,float roughness,float anistropic,float sheen,float sheenTint,float clearcoat,float clearcoartGlss){
    vec3 n = vec3(0.0,1.0,0.0);
    wi = sampleHemisphere(random.x,random.y,pdf);
    vec3 L = wo;
    vec3 V = wi;
    roughness = max(roughness * roughness,0.001);
    float NdotL = max(dot(n,wi),0.0);
    float NdotV = max(dot(n,wo),0.0);
    vec3 H = normalize(L + V);
    float NdotH = max(dot(n,H),0.0);
    float LdotH = dot(L,H); 

    //rho
    float luminance = 0.3 * baseColor.x + 0.6 * baseColor.y + 0.1 * baseColor.z;
    vec3 rho_tint = (luminance > 0.0) ? baseColor / luminance : vec3(1.0);
    vec3 rho_spec = lerp(vec3(1.0),rho_tint,specularTint);
    vec3 rho_sheen = lerp(vec3(1.0),rho_tint,sheenTint);

    //Diffse fresnel
    float Fd90 = 0.5 + 2.0 * LdotH* LdotH * roughness;
    float FL = SchlickFresnel(NdotL,Fd90),FV = SchlickFresnel(NdotV,Fd90);
    vec3 FD = baseColor * FL * FV / M_PI;

    //specular 
    float Fss90 = roughness * LdotH * LdotH;
    vec3 Fs0 = lerp(0.08 * specular * rho_spec,baseColor,metallic);
    vec3 mf = fresnelSchlick(Fs0,LdotH);
    float md = DGGX(NdotH,roughness); 
    float mg = GSJMSF(NdotV,NdotL,roughness);
    vec3 Fspec = Spec(wo,wi,roughness,Fs0);

    //sheen
    vec3 Fsheen = sheen * rho_sheen * pow(1.0 -LdotH,5.0);

    //subsurface
    float sFL = SchlickFresnel(NdotL,Fss90),sFV = SchlickFresnel(NdotV,Fss90);
    float subcos = 1. / (NdotL + NdotV) - 0.5;
    vec3 Fsub = baseColor * 1.25 * ( sFL * sFV * subcos + 0.5) /M_PI;

    //Clearcoat
    float alphac = lerp(0.1,0.001,clearcoartGlss);
    vec3 Fclear = clearC(wo,wi,alphac,clearcoat);

    return (lerp(FD,Fsub,subsurface)+ Fsheen)* (1.0 - metallic) + Fspec + Fclear;
}

//パストレーサー本体
vec3 pathtracer(Ray ray,vec2 uv){
	float p = 0.99;
	int Maxref = 10;
	vec3 LTE = vec3(0.0);
	vec3 s = vec3(1.0);
	Ray r = ray;

	float counter = 0.;
    vec3 direction = vec3(0.);
	for(int i = 0; i < 10; i++){

		if(random(uv * iTime) > p){
			break;
		}

		s /= p;

		Info info;
        info.distance = 1000.;

		if(!Scene(r,info)){
            //vec3 l = textureCube(cubeMap, vec3, 1.0)
			LTE = s * WorldSphere(r);
			break;
		}

		vec3 t,b;
        tangentSpaceBasis(info.normal,t,b);
		vec3 wo = worldtoLoacal(-ray.direction,t,info.normal,b);
		float pdf = 0.0;
        vec3 hash = vec3(hash22(uv),hash22(vec2(float(iTime) * 0.3,iTime)).x) + float(i) * 500.0 + 50.0;
        vec2 xi = hash23(hash);
        vec3 wi = vec3(1.0);
        
        // Define BSDF of each object
        // 0: floor
        // 1~5: index of each sphere from left;
        vec3 BSDF = vec3(0.0);
        if(info.number == 0){
            BSDF = Lambert(vec3(1.0),wi,pdf,xi);
        }
        else if(info.number  == 1){
            //BSDF = Lambert(vec3(1.0,0.0,0.0),wi,pdf,xi);
            BSDF = SmoothGlass(wo,wi,pdf,xi.x);
		}
        else if(info.number == 2){
             BSDF = Lambert(texture2D(texture02,info.uv).xyz,wi,pdf,xi);
        }
        else if(info.number == 3){
            BSDF = Lambert(texture2D(texture03,info.uv).xyz,wi,pdf,xi);
            //LTE = s * vec3(1.2,0.8,0.8) *2.0;
            //break;
        }
        else if(info.number == 4){
            BSDF = Lambert(texture2D(texture04,info.uv).xyz,wi,pdf,xi);
        }
        else if(info.number == 5){
            //BSDF = Lambert(vec3(0.0,0.0,1.0),wi,pdf,xi);
            
            BSDF = Mirror(vec3(0.9),wo,wi,pdf);
        }
        else if(info.number == 6){
            vec3 baseColor = _BASECOLOR/255.0;
            float subsurface = _SUBSURFACE;
            float metallic = _METALLIC;
            float specular = _SPECULAR;
            float specularTint = _SPECULARTINT;
            float roughness = _ROUGHNESS;
            float anistropic = _ANISTROPIC;
            float sheen = _SHEEN;
            float sheenTint = _SHEENTINT;
            float clearcoat = _CLEARCOAT;
            float clearcoartGlss = _CLEARCOAT;
            //BSDF = DisneyBRDF(wo, wi,xi,pdf,baseColor,subsurface, metallic,specular,specularTint,roughness,anistropic,sheen,sheenTint,clearcoat,clearcoartGlss);
            LTE = s * vec3(1.2,0.8,0.8) * _LUMINANCE;
            break;
        }
        else if(info.number == 10){
            //BSDF = Lambert(texture2D(texture01,info.uv).xyz,wi,pdf,xi);
            BSDF = MicroFacet(vec3(1.,0.86,0.57),wo,wi,clamp(texture2D(texture04,info.uv).x + 0.1,0.0,1.0),pdf,xi);
            //BSDF = MicroFacet(vec3(1.,0.86,0.57),wo,wi,0.01,pdf,xi);
        }
        else{
            BSDF = Lambert(texture2D(texture03,info.uv).xyz,wi,pdf,xi);
        }

		vec3 nextDirection = localToWorld(wi,t,info.normal,b);
        direction = nextDirection;
		float cosine = abs(dot(info.normal,nextDirection));
        s = s * BSDF * cosine / pdf;
		r = Ray(info.position + nextDirection * 0.01,nextDirection);
	}

	return LTE;

}

//カメラ基底を求める
mat3 camera(vec3 ro, vec3 ta)
{
    vec3 up = vec3(0, 1, 0);
    vec3 cw = normalize(ta - ro);
    vec3 cu = normalize(cross(cw, up));
    vec3 cv = normalize(cross(cu, cw));
    return mat3(cu, cv, cw);
}

mat3 cameraMat(vec3 ro)
{
    vec3 up = vec3(0, 1, 0);
    vec3 cw = normalize(ro);
    vec3 cu = normalize(cross(cw, up));
    vec3 cv = normalize(cross(cu, cw));
    return mat3(cu, cv, cw);
}

Ray thinLensCamera(vec2 uv,vec3 atlook,vec3 camerapos){
    float f = cameraLens.x;
    float F = cameraLens.y;
    float L = cameraLens.z;

    float V = L * f / (L - f);
    
    vec3 up = vec3(0,1,0);
    vec3 cw = normalize(atlook);
    vec3 cu = normalize(cross(cw,up));
    vec3 cv = normalize(cross(cu,cw));
    
    vec3 X = camerapos + uv.x * cu + uv.y * cv;
    vec3 C = camerapos + cw * V;
    vec3 e = normalize(C - X);

    //レンズ上の位置
    //ランダム関数
    vec2 xi = hash23(vec3(iTime * uv.x, iTime * iTime * uv.y , iTime * iTime * iTime));
    xi = hash22(xi);
    float phi = 2.0 * M_PI * xi.x;
    float r = xi.y * f / (2.0 * F);
    vec3 S = C + r * cos(phi) * cu + r * sin(phi) * cv;

    vec3 P = C + e * L / dot(e,cw);

    Ray camera;
    camera.origin = S;
    camera.direction = normalize(P - S);

    return camera;
}

vec3 normaltest(Ray r){
    Info info;
    Scene(r,info);
    return info.normal;
}
//出力部分
void main(void){
    vec2 hash = hash23(vec3(iTime * gl_FragColor.xy,iTime));
    float rand1 = random(hash) * 2.0 -1.0;
    float rand2 = random(hash.yx) * 2.0 -1.0;
    vec2 uv = (( gl_FragCoord.xy + vec2(rand1,rand2)) * 2.0 - iResolution.xy) / iResolution.y;
    vec2 lensuv = (((iResolution.xy - gl_FragCoord.xy) + vec2(rand1,rand2)) * 2.0 - iResolution.xy) / iResolution.y;
    Ray cameraray;    
    vec2 uv1 = gl_FragCoord.xy/iResolution.xy;
    
    //vec4 rayTip = cameraWorldMatrix * cameraProjectionInverseMatrix * vec4(uv.xy,1.,1.);
    vec3 rayTip = cameraMat(cameraDir) * vec3(uv,2.0);
    if(LensCheck){
        cameraray =  thinLensCamera(lensuv,cameraDir,cameraPos);
    }
    else{
        cameraray.origin = cameraPos;
        cameraray.direction = normalize(rayTip.xyz);
    }
    vec3 col = normaltest(cameraray);
    //今回の輝度を計算する
    if(!Normaltest){
        col =  pathtracer(cameraray,gl_FragCoord.xy); 
    }else{
        col = normaltest(cameraray);
    }
    vec4 prevColor = texture2D(buffer, uv1);
    if(iTime > 1.0){
        col += prevColor.xyz;
    }
    //前の出力に今回の輝度を加え、累計の輝度を出力す。。
    gl_FragColor = vec4(col,1.0);
}