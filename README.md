# deoptim

A Clojure library designed to ... well, that part is up to you.

## Usage

```clojure
(use 'deoptim.core)

(defn sphere [^doubles xs]
    (let [^doubles squares (amap ^doubles xs idx ret
                                 (Math/pow (aget ^doubles xs idx) 2))]
      (areduce ^doubles squares idx ret (double 0)
               (+ ret (aget ^doubles squares idx)))))

(optimize sphere 0.2 0.1 100 10 -100 100 150)
```

## License

Copyright © 2013 FIXME

Distributed under the Eclipse Public License, the same as Clojure.
