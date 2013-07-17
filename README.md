# deoptim

[![Build Status](https://secure.travis-ci.org/sorenmacbeth/deoptim.png?branch=master)](http://travis-ci.org/sorenmacbeth/deoptim)

A Clojure implementation of the Differential Evolution global
optimization algorithm.

To use in your leiningen projects add:

```clojure
[deoptim "0.1.0"]
```

to your project.clj

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

Copyright © 2013 Soren Macbeth

Distributed under the Eclipse Public License, the same as Clojure.
