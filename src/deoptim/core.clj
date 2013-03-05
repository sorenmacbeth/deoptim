(ns deoptim.core
  (:import [java.util Random]
           [org.apache.commons.math3.distribution NormalDistribution CauchyDistribution]))

;; ## An implementation of JADE
;; See http://sci2s.ugr.es/EAMHCO/pdf/JADE.pdf for details


;; # A record to represent a member of the population.

(defrecord Candidate [values crossover-prob mutation-factor fitness])

(def RANDOM
  (Random.))

(defn random-double
  [^Random r]
  (.nextDouble r))

(defn random-int
  [^Random r x]
  (.nextInt r x))

(defn next-double
  "Generate a double between start and end drawn from
   a uniform distribution."
  ^double [^double min ^double max]
  (let [^double d (random-double RANDOM)]
    (+ min (* d (- max min)))))

(defn doubles-seq
  "Infinite lazy-seq of doubles between min and max from
   a uniform distribution."
  [^double min ^double max]
  (lazy-seq
   (when true
     (cons (next-double min max) (doubles-seq min max)))))

;; # double array operations

(defn darr*
  "Multiply each element of `arr` by `f`."
  [^doubles arr ^double f]
  (amap ^doubles arr idx ret
        (* f (aget ^doubles arr idx))))

(defn darr-
  "Element-wise Subtraction of double array `y`
   from double array `x`."
  [^doubles x ^doubles y]
  (amap ^doubles x idx ret
        (- (aget ^doubles x idx) (aget ^doubles y idx))))

(defn darr+
  "Element-wise addition of double array `y`
   to double array `x`."
  [^doubles x ^doubles y]
  (amap ^doubles x idx ret
        (+ (aget ^doubles x idx) (aget ^doubles y idx))))

(defn darr=
  "Check if double arrays are equal."
  ([x] true)
  ([^doubles x ^doubles y] (java.util.Arrays/equals x y))
  ([x y & more]
     (if (darr= x y)
       (if (next more)
         (recur y (first more) (next more))
         (darr= y (first more)))
       false)))

(defn init-values
  "Create a single doubles vector of dimension `d` with values
   randomly chosen between `lower-bound` and `upper-bound`."
  ^doubles [d lower-bound upper-bound]
  (double-array (take d (doubles-seq lower-bound upper-bound))))

(defn init-population
  "Generate initial random population of size `size`
   as a seq of Candidate records."
  [size d ^double lower-bound ^double upper-bound]
  (into [] (for [_ (range size)
                 :let [^doubles v (init-values d lower-bound upper-bound)
                       cr 0.5
                       f 0.5]]
             (Candidate. v cr f Double/POSITIVE_INFINITY))))

(defn choose-p-best
  "Choose a random Candidate c from the p-best candidates
   in generation."
  [p c generation]
  (let [np (count generation)
        p-bests (int (Math/floor (* (* 1 p) np)))
        p-best (nth generation (random-int RANDOM p-bests))]
    (loop [p-best p-best]
      (if-not (darr= (:values c) (:values p-best))
        p-best
        (recur (nth generation (random-int RANDOM p-bests)))))))

(defn choose-random
  "Choose a random Candidate from generation."
  [c generation]
  (let [np (count generation)
        r (nth generation (random-int RANDOM np))]
    (loop [r r]
      (if-not (darr= (:values c) (:values r))
        r
        (recur (nth generation (random-int RANDOM np)))))))

(defn lehmer-mean
  "Calculate the Lehmer mean of `coll`."
  [p coll]
  (/ (reduce + (map #(Math/pow % p) coll))
     (reduce + (map #(Math/pow % (dec p)) coll))))

(defn update-mu-cr
  "Update `mu-cr` using constant `c`."
  [c mu-cr cr-set]
  (let [mean-cr (lehmer-mean 1 cr-set)]
    (+ (* (- 1 c) mu-cr) (* c mean-cr))))

(defn update-mu-f
  "Update `mu-f` using constant `c`."
  [c mu-f f-set]
  (let [mean-f (lehmer-mean 2 f-set)]
    (+ (* (- 1 c) mu-f) (* c mean-f))))

(defn update-cr
  "Update the crossover probability of a Candidate."
  [^NormalDistribution ndist c]
  (assoc c :crossover-prob (.sample ndist)))

(defn update-f
  "Update the mutation factor of a Candidate."
  [^CauchyDistribution cdist c]
  (assoc c :mutation-factor (loop [sample (.sample cdist)]
                              (cond
                               (> sample 1) 1.0
                               (<= sample 0) (recur (.sample cdist))
                               :else sample))))

(defn mutate
  "Generate a mutated double array based off of
   current values."
  [p generation archive c]
  (let [^double f (:mutation-factor c)
        ^doubles x (:values c)
        ^doubles p-best (:values (choose-p-best p c generation))
        ^doubles x1 (:values (choose-random c generation))
        ^doubles x2 (:values (choose-random c generation))]
    (assoc c :trial-values (darr+ x (darr+ (darr* (darr- p-best x) f) (darr* (darr- x1 x2) f))))))

(defn crossover
  "Perform crossover of mutated values.
   Ensure that `lower-bound` and `upper-bound`
   are respected."
  [jrand lower-bound upper-bound c]
  (let [^doubles x (:values c)
        ^doubles v (:trial-values c)
        crossover-prob (:crossover-prob c)]
    (assoc c :trial-values
           (into-array Double/TYPE
                       (map
                        (fn [idx value trial-value]
                          (if (or (< (random-double RANDOM) crossover-prob)
                                  (= idx jrand))
                            (cond
                             (< trial-value lower-bound) (/ (+ value lower-bound) 2)
                             (> trial-value upper-bound) (/ (+ value upper-bound) 2)
                             :else trial-value)
                            value))
                        (range (count x)) ; normal map-indexed only takes 1 coll arg
                        x v)))))

(defn selection
  "Perform selection amongst candidate and trial arrays."
  [fitness-fn archive c]
  (let [^double fit (fitness-fn (:values c))
        ^double trial-fit (fitness-fn (:trial-values c))]
    (if (< trial-fit fit)
      (-> (assoc c :values (:trial-values c) :fitness trial-fit) (dissoc :trial-values))
      (-> (assoc c :fitness fit) (dissoc :trial-values)))))

(defn evolve-candidate*
  "Evolve a single Candidate from a generation."
  [fitness-fn generation archive ndist cdist p jrand lower-bound upper-bound c]
  (->> (update-cr ndist c)
       (update-f cdist)
       (mutate p generation archive)
       (crossover jrand lower-bound upper-bound)
       (selection fitness-fn archive)))

(defn evolve-candidate [fitness-fn generation archive ndist cdist p jrand lower-bound upper-bound]
  (partial evolve-candidate* fitness-fn generation archive ndist cdist p jrand lower-bound upper-bound))

;; fitness-fn should accept a single double array as input and
;; return a double.

;; TODO: implement archive
(defn optimize
  "Do optimize the supplied function `fitness-fn` using the
   differential evolution algorithm."
  [fitness-fn p c np d lower-bound upper-bound maxiter]
  (loop [generation (init-population np d lower-bound upper-bound)
         mu-cr 0.5
         mu-f 0.5
         archive []
         maxiter maxiter]
    (let [^NormalDistribution ndist (NormalDistribution. mu-cr 0.1)
          ^CauchyDistribution cdist (CauchyDistribution. mu-f 0.1)
          jrand (random-int RANDOM d)
          new-generation (into [] (pmap
                                   (evolve-candidate fitness-fn generation archive ndist cdist p jrand lower-bound upper-bound) generation))
          sorted-generation (sort-by :fitness new-generation)]
      (if (neg? maxiter)
        (first sorted-generation)
        (do
          #_(println (str "[" maxiter "] best fitness: " (:fitness (first sorted-generation))))
          (recur sorted-generation
                (update-mu-cr c mu-cr (into #{} (map :crossover-prob new-generation)))
                (update-mu-f c mu-f (into #{} (map :mutation-factor new-generation)))
                []
                (dec maxiter)))))))

(comment
  (defn sphere [^doubles xs]
    (let [^doubles squares (amap ^doubles xs idx ret
                                 (Math/pow (aget ^doubles xs idx) 2))]
      (areduce ^doubles squares idx ret (double 0)
               (+ ret (aget ^doubles squares idx)))))

  (optimize sphere 0.2 0.1 100 10 -100 100 150))
