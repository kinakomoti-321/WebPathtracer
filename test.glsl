#define M_PI float(3.141592)
#define HASHSCALE3 vec3(.1031, .1030, .0973)
#define lerp(v1,v2,x) x * v1 + (1.0 - x) * v2
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

//ランダム関数、サンプリングなどに使う。
float random (vec2 st) {
    return fract(dot(st.xy,
                         vec2(12.9898,78.233)));
}

//reffer https://github.com/gam0022/webgl/blob/master/pathtracing_sandbox.html by gam0022
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
    if(t < 0.0001 || t > 10000.){
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
bool baseplaneIntersect(Ray ray,inout Info info){
    vec3 normal = normalize(vec3(0.,1.,0.));
    
    float sn = dot(ray.origin,normal);
    float dn = dot(ray.direction,normal);

    if(abs(dn) == 0.000) return false;

    float t = - sn/ dn;
    if(t < 0.0001 || t > 10000.) return false;
    vec3 pos = ray.origin + t * ray.direction;
    if(length(pos) > 8.) return false;
    info.distance = t;
    info.normal = normal;
    info.position = pos;
    info.color = vec3(0.9);
    info.number = 0;
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
    vec3 radiance =  texture2D(texture01,worlduv).xyz;
    //radiance = vec3(0.0);
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
    if(baseplaneIntersect(r,tes)){
        if(tes.distance < info.distance){
            hit = true;
            info = tes;
        }
    }
    return hit;
}

//BRDF
vec3 Lambert(vec3 rho,inout vec3 wi,inout float pdf,vec2 random){
     wi = sampleHemisphere(random.x,random.y,pdf);
     return rho / M_PI;
}
//MicrofasetBRDF
// https://qiita.com/aa_debdeb/items/f813bdcbd8524a66a11b
vec3 fresnelSchlick(vec3 f0,float cosine){
    return f0 + (1.0 - f0) * pow(1.0 - cosine,5.0);
}

float DGGX(float nh,float alpha){
    float d = (1.0 - (1.0 - alpha) * nh * nh);
    return alpha/(M_PI * d * d);
}

float Lambda(float alpha,float xn){
    return 0.5 * (-1.0 + sqrt(1.0 + alpha * (1./(xn * xn) - 1.0)));
}

float GSJMSF(float nv,float nl,float alpha){
    float lambdav = Lambda(alpha,nv);
    float lambdal = Lambda(alpha,nl);
    return 1./(1. + lambdav + lambdal);
}


vec3 MicroFacet(vec3 F0, vec3 wo,inout vec3 wi,float roughness,inout float pdf,vec2 random){
    vec3 n = vec3(0.0,1.0,0.0);
    wi = sampleHemisphere(random.x,random.y,pdf);
    vec3 v = wo;
    vec3 l = wi;
    vec3 h = normalize(v+l);
    float nh = dot(n,h);
    float nl = dot(n,l);
    float nv = dot(n,v);
    float vh = dot(v,h);
    
    float alpha = roughness * roughness;
    
    float D = DGGX(nh,alpha);
    float G = GSJMSF(nv,nl,alpha);
    vec3 F = fresnelSchlick(F0,vh);
    
    vec3 f = F * D * G / (4.0 * nv * nl);
    
    return f;

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
        vec3 hash = vec3(uv,float(iTime) * 0.3) + float(i) * 500.0 + 50.0;
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
            BSDF = Lambert(vec3(1.0,0.0,0.0),wi,pdf,xi);
		}
        else if(info.number == 2){
             BSDF = Lambert(texture2D(texture02,info.uv).xyz,wi,pdf,xi);
        }
        else if(info.number == 3){
            // BSDF = Lambert(texture2D(texture03,info.uv).xyz,wi,pdf,xi);
            LTE = s * vec3(1.2,0.8,0.8) *2.0;
            break;
        }
        else if(info.number == 4){
            BSDF = Lambert(texture2D(texture03,info.uv).xyz,wi,pdf,xi);
        }
        else if(info.number == 5){
            BSDF = Lambert(vec3(0.0,0.0,1.0),wi,pdf,xi);
        }
        else{
            BSDF = vec3(0.0);
            wi = vec3(0.0,1.0,0.0);
            pdf = 1.;
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
    float f = 2.0;
    float F = 2.0;
    float L = 1.0;

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
    float phi = 2.0 * M_PI * xi.x;
    float r = xi.y * f / (2.0 * F);
    vec3 S = C + r * cos(phi) * cu + r * sin(phi) * cv;

    vec3 P = C + e * L / dot(e,cw);

    Ray camera;
    camera.origin = S;
    camera.direction = normalize(P - S);

    return camera;
}

//出力部分
void main(void){
    vec2 hash = hash23(vec3(iTime * gl_FragColor.xy,iTime));
    float rand1 = random(hash) * 2.0 -1.0;
    float rand2 = random(hash.yx) * 2.0 -1.0;
    vec2 uv = ((gl_FragCoord.xy + vec2(rand1,rand2)) * 2.0 - iResolution.xy) / iResolution.y;
    Ray cameraray;    
    vec2 uv1 = gl_FragCoord.xy/iResolution.xy;
    
    //vec4 rayTip = cameraWorldMatrix * cameraProjectionInverseMatrix * vec4(uv.xy,1.,1.);
    vec3 rayTip = cameraMat(cameraDir) * vec3(uv,2.0);
    cameraray.origin = cameraPos;
    cameraray.direction = normalize(rayTip.xyz);
    //今回の輝度を計算する
    vec3 col =  pathtracer(cameraray,gl_FragCoord.xy); 
    vec4 prevColor = texture2D(buffer, uv1);
    if(iTime > 1.0){
        col += prevColor.xyz;
    }
    //前の出力に今回の輝度を加え、累計の輝度を出力す。。
    gl_FragColor = vec4(col,1.0);
}