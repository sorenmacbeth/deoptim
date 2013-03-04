(ns deoptim.core
  (:import [java.util Random]
           [org.apache.commons.math3.distribution NormalDistribution CauchyDistribution]))

;; ## An implementation of JADE
;; See http://sci2s.ugr.es/EAMHCO/pdf/JADE.pdf for details


;; # A record to represent a member of the population.

(defrecord Candidate [values crossover-prob mutation-factor fitness])

(def RANDOM
  (Random.))

(defn next-double
  "Generate a double between start and end drawn from
   a uniform distribution."
  ^double [^double min ^double max]
  (let [^double d (.nextDouble RANDOM)]
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
        p-best (nth generation (.nextInt RANDOM p-bests))]
    (loop [p-best p-best]
      (if-not (darr= (:values c) (:values p-best))
        p-best
        (recur (nth generation (.nextInt RANDOM p-bests)))))))

(defn choose-random
  "Choose a random Candidate from generation."
  [c generation]
  (let [np (count generation)
        r (nth generation (.nextInt RANDOM np))]
    (loop [r r]
      (if-not (darr= (:values c) (:values r))
        r
        (recur (nth generation (.nextInt RANDOM np)))))))

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

(defn generate-trial-values
  "Generate a double array of trial values."
  [p generation c]
  (let [^double f (:mutation-factor c)
        ^doubles x (:values c)
        ^doubles p-best (:values (choose-p-best p c generation))
        ^doubles x1 (:values (choose-random c generation))
        ^doubles x2 (:values (choose-random c generation))]
    (assoc c :trial-values (darr+ x (darr+ (darr* (darr- p-best x) f) (darr* (darr- x1 x2) f))))))

(defn mutate
  "Mutate a candidate record."
  [jrand c]
  (let [^doubles u (:values c)
        ^doubles v (:trial-values c)]
    (assoc c :values (amap ^doubles u idx ret
                           (if (or (= idx jrand) (< (.nextDouble RANDOM) (:crossover-prob c)))
                             (aget ^doubles v idx)
                             (aget ^doubles u idx))))))

(defn selection
  "Perform selection amongst candidate and trial arrays."
  [fitness-fn c]
  (let [fit (fitness-fn (:values c))
        trial-fit (fitness-fn (:trial-values c))]
    (if (< trial-fit fit)
      (-> (assoc c :values (:trial-values c) :fitness trial-fit) (dissoc :trial-values))
      (-> (assoc c :fitness fit) (dissoc :trial-values)))))

;; fitness-fn should accept a single double array as input and
;; return a double.
(defn optimize
  "Do optimize the supplied function `fitness-fn` using the
   differential evolution algorithm."
  [fitness-fn p c np d lower-bound upper-bound maxiter]
  (loop [generation (init-population np d lower-bound upper-bound)
         mu-cr 0.5
         mu-f 0.5
         archive []
         maxiter maxiter]
    (let [norm (NormalDistribution. mu-cr 0.1)
          cauchy (CauchyDistribution. mu-f 0.1)
          jrand (.nextInt RANDOM d)
          new-generation (->>
                          (pmap #(update-in % [:crossover-prob]
                                            (fn [_]
                                              (.sample norm))) generation)
                          (pmap #(update-in % [:mutation-factor]
                                            (fn [_]
                                              (loop [sample (.sample cauchy)]
                                                (cond
                                                 (> sample 1) 1.0
                                                 (<= sample 0) (recur (.sample cauchy))
                                                 :else sample)))))
                          (pmap #(generate-trial-values p generation %))
                          (pmap #(mutate jrand %))
                          (pmap #(selection fitness-fn %)))]
      (if (zero? maxiter)
        (first (sort-by :fitness #(compare %2 %1) new-generation))
        (recur (sort-by :fitness #(compare %2 %1) new-generation)
               (update-mu-cr c mu-cr (into #{} (map :crossover-prob new-generation)))
               (update-mu-f c mu-f (into #{} (map :mutation-factor new-generation)))
               []
               (dec maxiter))))))