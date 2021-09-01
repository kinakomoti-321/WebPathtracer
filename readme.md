# web RayTracer
Three.jsを使用したウェブ上で動くレイトレーサーです。

PC、Chromeの使用を推奨します。（現在画面が４分割されるようなバグが確認されています。環境によってはこのバグが生じる可能性があります。）

DemoPage 
https://kinakomoti-321.github.io/WebPathtracer/

## GUIの各パラメーター
1. Camera

カメラに関するGUIです。

	- f : カメラの焦点距離
	- F : カメラのF値　大きいほどピントが外れてもぼやけにくい
	- L : レンズからピントまでの距離
	- Lens : 薄レンズのON,OFF
	- Shot : スクリーンショット

2. Scene

シーンに関するGUIです。

	- Background : 背景のON,OFF

3. Material

Disney BRDFを適用した球体に関するGUIです。

	- BaseColor : 物体の色
	- Subsurface :　表面化散乱の度合い
	- Roughness :　粗さの度合
	- Specular :　非金属のスペキュラー反射の調節
	- SpecularTint :　スペキュラー反射の色がBaseColorを基にするかの度合
	- Matellic : 金属の度合い
	- Anistropic :　異方性
	- Sheen :　布っぽくするパラメータ
	- SheenTint :　Sheenによる反射光がBaseColorの色になるかの度合
	- Clearcoat :　滑らかなニスを塗ったような反射を付ける項目
	- ClearcoatTint :　Clearcoatの反射がBaseColorの色になるかの度合い
	- x : x座標
	- y : y座標
	- z : z座標

画像の取得先

earth.jpg
[NASA visible earth](https://visibleearth.nasa.gov/images/57730/the-blue-marble-land-surface-ocean-color-and-sea-ice/57731l)

sun.png [https://www.textures.com/download/PBR0221/133269]

参考文献
@gamm0022
[WebGL+GLSLによる超高速なパストレーシング](https://qiita.com/gam0022/items/18bb3612d7bdb6f4360a)

Pent@creation Blog [Three.jsでオフスクリーンレンダリング](https://www.pentacreation.com/blog/2021/05/210504.html)

@kaneta1992 
[three.js + キューブマップでお手軽IBL](https://qiita.com/kaneta1992/items/df1ae53e352f6813e0cd)

