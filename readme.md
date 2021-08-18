# web RayTracer
Three.jsを使用したウェブ上で動くレイトレーサーです。

PC、Chromeの使用を推奨します。（現在画面が４分割されるようなバグが確認されています。環境によってはこのバグが生じる可能性があります。）

$ \frac{G F}{4|\omega_i \cdot n| |\omega_m \cdot n|} $

DemoPage 
https://kinakomoti-321.github.io/WebPathtracer/

## GUIの各パラメーター
1. Camera
	- f : カメラの焦点距離
	- F : カメラのF値　大きいほどピントが外れてもぼやけにくい
	- L : レンズからピントまでの距離
	- Lens : 薄レンズのON,OFF
	- Shot : スクリーンショット

2. Scene
	- Background : 背景のON,OFF


参考文献

@gamm0022
[WebGL+GLSLによる超高速なパストレーシング](https://qiita.com/gam0022/items/18bb3612d7bdb6f4360a)

Pent@creation Blog [Three.jsでオフスクリーンレンダリング](https://www.pentacreation.com/blog/2021/05/210504.html)

@kaneta1992 
[three.js + キューブマップでお手軽IBL](https://qiita.com/kaneta1992/items/df1ae53e352f6813e0cd)

