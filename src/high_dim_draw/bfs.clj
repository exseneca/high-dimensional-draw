(ns high-dim-draw.bfs
  (:require [clojure.set :as set]
            [clojure.reflect :as r]
            [clojure.math :as m]
            [clatrix.core :as cl :exclude [abs]]
            [clojure2d.core :as draw]
            [clojure.core.matrix :as matrix]
            [clojure.core.matrix.linear :as linear]
            [clojure.pprint :as pp])
  (:import (java.awt.event KeyListener KeyEvent)))

;; graph is a collection of nodes and edges
(def graph [{:a [:b :h :d]}
            {:b [:a :c :d]}
            {:c [:b :e :d]}
            {:d [:b :a :c]}
            {:e [:c :f :g :h]}
            {:f [:e :g]}
            {:g [:e :f]}
            {:h [:e :a]}])

;(def big-graph-nodes (repeatedly 100 (fn [] (keyword (str (random-uuid))))))

;
;; (def graph (vec (->>  (into {} (map
;;                                 (fn [x]
;;                                   [x (vec (random-sample 0.1 big-graph-nodes))])
;;                                 big-graph-nodes))
;;                       (map (fn [[k v]] (hash-map k v))))))
(defn make-queue []
  (clojure.lang.PersistentQueue/EMPTY))

;; Parent needs to be added. Some more thought needs to be done on this
;; (defn bfs- [graph end-node visited to-visit]
;;   (if (empty? to-visit)
;;     (throw (ex-info "Can't find it!" {}))
;;     (let [next-node (peek to-visit)]
;;       (if (visited next-node)
;;         (bfs- graph end-node visited (pop to-visit))
;;         (let [new-paths next-node])))))

(defn get-next-node [queue]
  (peek queue))

(defn get-next-node-value [queue]
  (first (get-next-node queue)))

(defn get-node-value [node]
  (first node))

(defn get-node-distance [node]
  (second node))

(defn get-node-paths
  "From a graph which is a vector of maps.
  The map being a key of the node name and
  a value which is a vector of adjacent nodes,
  get the adjacent nodes for a given node."
  [graph node]
  (-> (get-node-value node)
      (filter graph)
      first
      (get (get-node-value node))))


(defn get-seen-paths-set [visited to-visit]
  (set/union (set (map get-node-value to-visit)) ;; visited already a set
             visited))
(declare bfs-)

(defn add-new-nodes [graph end-node visited to-visit node]
  (let [new-paths  (get-node-paths graph node)
        node-distance (get-node-distance node)]
    (if (some #(= % end-node) new-paths)
      (inc node-distance)
      (let [seen-paths (get-seen-paths-set visited to-visit)
            unseen-path-values (set/difference (set new-paths) seen-paths)
            unseen-nodes (map (fn [x y] [x y]) unseen-path-values (repeat (inc node-distance)))]
        (bfs- graph end-node (conj visited (get-node-value node)) (apply (partial conj (pop to-visit)) unseen-nodes))))))

(defn bfs- [graph end-node visited to-visit]
  (if (empty? to-visit)
    (throw (ex-info "Can't find it!" {}))
    (add-new-nodes graph end-node visited to-visit (peek to-visit))))

(defn bfs [graph start-node end-node]
  (if (= start-node end-node)
    0
    (bfs- graph end-node #{} (conj (make-queue) [start-node 0]))))

(defn test-bfs []
  (assert (= 2 (bfs graph :a :e))))

(defn key-from-node [node]
  (-> node
      first
      first))

(defn random-node [graph]
  (->> (count graph)
       (rand-int)
       (nth graph)
       (key-from-node)))

(defn get-next-pivot [d-vec nodes]
  (->> (apply max d-vec)
       (.indexOf d-vec)
       (nth nodes)))

(defn high-dim-draw- [graph nodes pivot out-vectors d-vec m]
  (if (= m 0)
    out-vectors
    ;; Could speed things up if we 
    (let [pivot-values (mapv (partial bfs graph pivot) nodes)
          new-d-vec (map min d-vec pivot-values)
          next-pivot (get-next-pivot new-d-vec nodes)]
      (high-dim-draw- graph nodes next-pivot (conj out-vectors pivot-values) new-d-vec (dec m)))))

;; This will break if the graph size is greater than Integer/MAX_VALUE
(defn high-dim-draw [graph dimension]
  (let [start-node (random-node graph)
        nodes (map key-from-node graph)
        d-vec (repeat Integer/MAX_VALUE (count graph))]
    (high-dim-draw- graph nodes start-node [] d-vec dimension)))

(defn covariance [m]
  (matrix/mmul m (matrix/transpose m)))

(defn eigenvectors [cov]
  (let [eigens (:Q (linear/eigen (cl/matrix cov)))]
    [(first eigens) (second eigens)]))

(defn eigenvectors-3d [cov]
  (let [eigens (:Q (linear/eigen (cl/matrix cov)))]
    [(first eigens) (second eigens) (nth eigens 2)]))

;; (def out-matrix (high-dim-draw graph 3))
;; ? (matrix/mmul (mapv vec (eigenvectors (covariance out-matrix))) out-matrix)

(defn to-vec ([x y]
              [x y])
  ([x y z] [x y z]))

(defn coords [graph dimension]
  (let [out-matrix (high-dim-draw graph dimension)]
    (println "eigenvectors: " (eigenvectors (covariance out-matrix)))
    (println "eigenvectors mapped: " (mapv vec (eigenvectors (covariance out-matrix))))
    (->> (matrix/mmul (mapv vec (eigenvectors (covariance out-matrix))) out-matrix)
         (apply (partial map to-vec)))))

(defn coords-3d [graph dimension]
  (let [out-matrix (high-dim-draw graph dimension)]
    (->> (matrix/mmul (mapv vec (eigenvectors-3d (covariance out-matrix))) out-matrix)
         (apply (partial map to-vec)))))

(defn scale-coords [coords]
  (let [xs (map first coords)
        ys (map second coords)
        x-min (apply min xs)
        x-max (apply max xs)
        y-min (apply min ys)
        y-max (apply max ys)
        x-diff (- x-max x-min)
        y-diff (- y-max y-min)
        scale (* 1 (+ 0.5 (max x-diff y-diff)))
        x-move (if (< x-min 0) (abs x-min) 0)
        y-move (if (< y-min 0) (abs y-min) 0)]
    (map (fn [[x y]] [(/ (+ x x-move 0.1) scale) (/ (+ y y-move 0.1) scale)]) coords)))

(defn scale-coords-3d [coords]
  (let [xs (map first coords)
        ys (map second coords)
        zs (map (fn [x] (nth x 2)) coords)
        x-min (apply min xs)
        x-max (apply max xs)
        y-min (apply min ys)
        y-max (apply max ys)
        z-max (apply max zs)
        z-min (apply min zs)
        x-diff (- x-max x-min)
        y-diff (- y-max y-min)
        z-diff (- z-max z-min)
        scale (* 1 (+ 0.5 (max z-diff x-diff y-diff)))
        ;; TODO: align this properly. 
        x-move (if (< x-min 0) (abs x-min) 0)
        y-move (if (< y-min 0) (abs y-min) 0)
        z-move (if (< z-min 0) (abs z-min) 0)]
    (map (fn [[x y z]] [(/ (+ x x-move) scale) (/ (+ y y-move) scale) (/ (+ z z-move) scale)]) coords)))

(def points (-> (coords graph 4)
                (scale-coords)))

(def points-3d
  (-> (coords-3d graph 4)
      (scale-coords-3d)))

(defn graph-keys [g] (map first (map first g)))
(defn point-lookup [pts] (apply hash-map (interleave (map first (map first graph)) pts)))

(defn expand-to-points [node-lookup [node connected-nodes]]
  (let [node-point (get node-lookup node)]
    (map (fn [x] [node-point (get node-lookup x)]) connected-nodes)))

(defn create-edges [graph points]
  (let [lookup (point-lookup points)]
    (->>
     (map  (partial expand-to-points lookup) (map first graph))
     (map first))))

(defn translate-to-canvas [x y]
  [(* 1000 x) (- 1000 (* 1000 y))])

(defn draw-line [canvas [[x1 y1] [x2 y2]]]
  (draw/set-color canvas :gray)
  (draw/line canvas (* 1000 x1) (- 1000 (* 1000 y1)) (* 1000 x2) (- 1000 (* 1000 y2))))
(defn draw-red-line [canvas [[x1 y1] [x2 y2]]]
  (draw/set-color canvas :red)
  (draw/line canvas (* 1000 x1) (- 1000 (* 1000 y1)) (* 1000 x2) (- 1000 (* 1000 y2))))

(defonce *window
  (atom nil))

(defonce canvas (draw/canvas 1000 1000))

(defn draw-fn-2d [points graph]
  (fn
    [canvas _ _ _]
    (-> canvas
        (draw/set-color 45 45 41 20))
    (doseq [point points]
      (draw/ellipse canvas (int (* 1000 (first point))) (- 1000 (int (* 1000 (second point)))) 10 10))
    (doseq [line (create-edges graph points)]
      (draw-line canvas line))))

;; (def window (draw/show-window {:canvas (draw/canvas 1000 1000)
;;                                :draw-fn (draw-fn points graph)}))

(defn scale-point [[x y] x-scale y-scale]
  [(* x-scale x) (* y-scale y)])

(defn translate-point [[x y] x-translate y-translate]
  [(+ x x-translate) (+ y y-translate)])

;; How much to translate the x 
(defn x-plane-translate [[base-line-x0 _] [base-line-x1 _] scale]
  (let [new-length (* scale (- base-line-x1 base-line-x0))]
    (/ (+ base-line-x1 base-line-x0 (- new-length)) 2)))

(defn translate-to-plane [point base-start base-end scale height]
  (-> point
      (scale-point scale scale)
      (translate-point (x-plane-translate base-start base-end scale)
                       (- height (* scale height)))))
(assert (= [25.0 50.0]
           (translate-to-plane [0 0] [0 0] [100 0] 0.5 100)))

(def base-line [[0.1 0.1] [0.9 0.1]])

(def next-line (mapv (fn [p]
                       (translate-to-plane p (first base-line) (second base-line) 0.5 1))
                     base-line))

(defn line-scale [base-line scale]
  (mapv (fn [p]
          (translate-to-plane p (first base-line) (second base-line) scale 0.8))
        base-line))
(def rotation-proportion 0.1)
(defn get-rotation-value [rotation-proportion] (* m/PI 2 rotation-proportion))

                                        ;                                   ;(def camera [0.5 0.5 -2])
;; range [2.5, infinity)
(defonce zoom-a (atom 0))

(defn camera-x [theta zoom]
  (+ 0.5 (* (- zoom 2.5)
            (m/sin theta))))
(defn camera-z [theta zoom]
  (+ 0.5 (* (- zoom 2.5)
            (m/cos theta))))

(defn get-camera [rotation-value zoom] [(camera-x rotation-value zoom) 0.5 (camera-z rotation-value zoom)])

(def camera-plane-distance 1.5)

(defn get-camera-orientation [rotation-value] [0 rotation-value 0])

(defn project-3d->2d [{:keys [camera camera-orientation]} point]
  (let [[theta-x theta-y theta-z] camera-orientation
        [dx dy dz]
        (->
         (matrix/mmul [[1 0 0] [0 (m/cos theta-x) (m/sin theta-x)]
                       [0 (- (m/sin theta-x)) (m/cos theta-x)]]
                      [[(m/cos theta-y) 0 (- (m/sin theta-y))]
                       [0 1 0]
                       [(m/sin theta-y) 0 (m/cos theta-y)]])
         (matrix/mmul [[(m/cos theta-z) (m/sin 0) theta-z]
                       [(- (m/sin theta-z)) (m/cos theta-z) 0]
                       [0 0 1]])
         (matrix/mmul (matrix/sub point camera)))]
    [(+ 0.5 (* camera-plane-distance (/ dx dz)))
     (+ 0.5 (* camera-plane-distance (/ dy dz)))]))

(defn translate [[x y z] d]
  [(Math/pow x (/ d (+ z d)))
   (Math/pow y (/ d (+ z d)))])

(defn add-z [points z]
  (mapv (fn [p] (conj p z)) points))

(defn project-line [ctx [point-start point-end :as line] z d]
  (mapv (fn [p]
          (project-3d->2d ctx p))
        (add-z line z)))

(defn project-line-3d [ctx points d]
  (mapv (fn [p] ;(translate p d)
          (project-3d->2d ctx p)) points))

(defonce rotation-a (atom 0))

(declare re-draw)
(def increment-amount 0.005)
(defn left-pressed []
  (swap! rotation-a #(- % increment-amount)))
(defn right-pressed []
  (swap! rotation-a #(+ % increment-amount)))

(defn plus-pressed []
  (swap! zoom-a (fn [zoom] (min 2.5 (+ zoom (* 10 increment-amount))))))
(defn minus-pressed []
  (swap! zoom-a (fn [zoom] (- zoom (* 10 increment-amount)))))

(defn make-key-pressed []
  (proxy [KeyListener] []
    (keyPressed [^KeyEvent key-press]
      (case (.getKeyCode key-press)
        37 (right-pressed)
        39 (left-pressed)
        61 (plus-pressed)
        45 (minus-pressed)
        nil))
    (keyReleased [^KeyEvent _])
    (keyTyped [^KeyEvent _])))

(defn make-context [rot zoom]
  (let [rotation-value (get-rotation-value rot)]
    {:camera (get-camera rotation-value zoom)
     :camera-orientation (get-camera-orientation rotation-value)}))

(defn draw-label [cv text x y]
  (draw/set-color cv java.awt.Color/BLACK)
  (let [[canvas-x canvas-y] (translate-to-canvas x y)]
    (draw/text cv (apply str (take 6 text))  canvas-x (+ (- 10) canvas-y) :left)))

;; TODO: refactor into a scene that could be translated
(defn draw-fn [rotation-a]
  (fn [canvas _ _ _]
    (let [ctx (make-context @rotation-a @zoom-a)]
      (-> canvas
          (draw/set-color :white)
          (draw/rect 0 0 1000 1000)
          (draw/set-color 45 45 41 20))
      (let [lines (map (fn [p] (project-line ctx base-line p 0.5)) [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9])]
;      (draw-line canvas base-line)
                                        ;      (last lines)

        (draw-line canvas (project-line-3d ctx [[0.1 0.1 0.1] [0.1 0.1 0.9]] 0.5))
        (doseq [h [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8]]
          (draw-line canvas (project-line-3d ctx [[(+ h 0.1) 0.1 0.1] [(+ h 0.1) 0.1 0.9]] 0.5)))
        (doseq [h [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8]]
          (draw-line canvas (project-line-3d ctx [[0.1 (+ h 0.1) 0.1] [0.1 (+ h 0.1) 0.9]] 0.5)))
        (doseq [h [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9]]
          (draw-line canvas (project-line-3d ctx [[0.1 0.1 h] [0.1 0.9 h]] 0.5)))
        (doseq [h [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8]]
          (draw-line canvas (project-line-3d ctx [[0.1 (+ h 0.1) 0.9] [0.9 (+ h 0.1) 0.9]] 0.5)))
        (doseq [h [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8]]
          (draw-line canvas (project-line-3d ctx [[(+ 0.1 h) 0.1 0.9] [(+ 0.1 h) 0.9 0.9]] 0.5)))
        (doseq [line  lines]
          (draw-line canvas line))
        (let [node-names (map (fn [x] (first (first x))) graph)
              name-and-points (map (fn [x y] [x y]) node-names points-3d)]
          (doseq [[label point] name-and-points]
            (let [translate-point (project-3d->2d ctx point)]
              (draw-label canvas (name label) (first translate-point) (second translate-point))

              (draw/ellipse canvas (int (* 1000 (first translate-point)))
                            (- 1000 (int (* 1000 (second translate-point)))) 10 10))))
        (doseq [line (create-edges graph points-3d)]
          (draw-red-line canvas (project-line-3d ctx line 0.5)))))))

(def window-2d (atom nil))

(def window (atom nil))

(defn draw []
  (reset! window-2d
          (draw/show-window {:canvas (draw/canvas 1000 1000)
                             :draw-fn (draw-fn-2d points graph)}))
  (reset! window
          (draw/show-window {:canvas (draw/canvas 1000 1000 :highest)
                             :draw-fn (draw-fn rotation-a)})))

(defn close []
  (.setVisible (:frame @window) false)
  (.dispose (:frame @window))
  (reset! window false))

(defn re-draw []
  (when @window (close))
  (draw)
  (.addKeyListener (:panel @window) (make-key-pressed)))


