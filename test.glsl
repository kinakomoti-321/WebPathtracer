
//Disney BRDF
vec3 disneyBRDF(){
    vec3 n = vec3(0.0,1.0,0.0);
    wi = sampleHemisphere(random.x,random.y,pdf);
    return vec3(0.0);
}