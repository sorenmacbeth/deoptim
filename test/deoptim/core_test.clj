(ns deoptim.core-test
  (:use clojure.test
        deoptim.core))

(deftest a-test
  (testing "Sphere is well optimized."
    (defn sphere [^doubles xs]
      (let [^doubles squares (amap ^doubles xs idx ret
                                   (Math/pow (aget ^doubles xs idx) 2))]
        (areduce ^doubles squares idx ret (double 0)
                 (+ ret (aget ^doubles squares idx)))))
    (time 
     (optimize sphere 0.2 0.1 100 10 -100 100 150))
    ;; initial trial times (ms): 5896, 5845, 5898, 5910, 6009
    ;; after rewriting crossover (ms): 2354, 2363, 2437, 2392, 2546    
    ;; after also using hinted random-* fns (ms): 1628, 1596, 1620, 1570, 1590
    ))
