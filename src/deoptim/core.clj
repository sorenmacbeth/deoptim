(ns deoptim.core
  (:require [bigml.sampling [simple :as simple]]))

;; ## An implementation of JADE
;; See http://sci2s.ugr.es/EAMHCO/pdf/JADE.pdf for details

(defn init-candidate
  "Create a single random candidate vector."
  ^doubles [d constraints]
  (double-array (take d (simple/sample constraints :replace true :generator :twister))))

(defn init-population
  "Generate initial random population of size `size`
   as a seq of double[]."
  [size d constraints]
  (for [_ (range size)
        :let [c (init-candidate d constraints)]]
    c))

;; TODO: Make this work.
(defn optimize
  "Do optimize the supplied function `eval-fn` using the
   differential evolution algorithm."
  [eval-fn & {:keys [np f cr stategy dither max-iter] :as opts}]
  (loop [population (init-population np) iter 0]
    (let [vals (map eval-fn population)
          done? false]
      (if (or done? (> iter max-iter))
        :return-best-candidate
        (recur population (inc iter))))))