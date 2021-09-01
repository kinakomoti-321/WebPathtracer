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

	- Background : 
		- World: 背景あり
		- NoBackground: 背景なし
		- WhiteFurnanceTest: 白背景
	- normaltest : ノーマルを可視化する

3. Light

ライトの球に関するGUIです。

	- x : x座標
	- y : y座標
	- z : z座標
	- Radius : 球の半径
	- Luminance : 光度

画像の取得先

earth.jpg
[NASA visible earth](https://visibleearth.nasa.gov/images/57730/the-blue-marble-land-surface-ocean-color-and-sea-ice/57731l)

Mars.jpg [NASA 3D resources](https://nasa3d.arc.nasa.gov/detail/mar0kuu2)

moon.jpg [NASA Scientific Visualization Studio](https://svs.gsfc.nasa.gov/4720)

参考文献
@gamm0022
[WebGL+GLSLによる超高速なパストレーシング](https://qiita.com/gam0022/items/18bb3612d7bdb6f4360a)

Pent@creation Blog [Three.jsでオフスクリーンレンダリング](https://www.pentacreation.com/blog/2021/05/210504.html)

@kaneta1992 
[three.js + キューブマップでお手軽IBL](https://qiita.com/kaneta1992/items/df1ae53e352f6813e0cd)

