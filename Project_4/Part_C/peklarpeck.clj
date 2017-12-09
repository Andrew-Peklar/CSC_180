 (ns fungp.peklarpeck
  (:use fungp.core)
  (:use fungp.util)
  (:use clojure.pprint))

(def criteria 0.01)

(def sample-functions
  '[[+ 2]
    [- 2]
    [* 2]
    [fungp.util/abs 1]
    [fungp.util/cos 1]
    [fungp.util/sdiv 2]
    [fungp.util/sin 1]
    [fungp.util/sqrt 1]
    [inc 1]
    [dec 1]])

(def sample-parameters
    ['x 'y])

(def number-literals
  (map float (range 10)))

(def in-list1 '(-25 -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25))
(def in-list2 '(-25 -24 -23 -22 -21 -20 -19 -18 -17 -16 -15 -14 -13 -12 -11 -10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25))
(def out-list '(50.52474970740786 49.53650932264733 44.19642330470238 43.96459614978917 43.833043095831265 38.5097736790413 37.40726284258123 37.98355770688623 32.94183462775995 30.89714663751662 31.976063248185724 27.45818842338426 24.474883099040795 25.811156724013248 22.017702618580806 18.174109498544745 19.50197449354335 16.57580663333013 12.018785288610259 13.07314583600087 11.088042221778739 6.021283506753236 6.558830996397852 5.513604990615857 0.18140514634863658 0.0 3.8185948536513634 2.4863950093841436 5.441169003602148 9.978716493246763 8.911957778221261 10.92685416399913 15.981214711389741 15.424193366669869 16.49802550645665 21.825890501455255 21.982297381419194 22.188843275986752 27.525116900959205 28.54181157661574 28.023936751814276 33.10285336248338 35.05816537224005 34.01644229311377 38.59273715741877 41.4902263209587 40.166956904168735 44.03540385021083 47.80357669529762 46.46349067735267 49.47525029259214))

(defn sample-fitness
  [tree]
  (try
    (let [f (compile-tree tree sample-parameters)
          results (map f in-list1 in-list2)]
      (reduce + (map off-by-sq out-list results)))
    (catch Exception e (println e) (println tree))))

(defn sample-report
  [tree fitness]
  (pprint tree)
  (println (str "Error:\t" fitness "\n"))
  (flush))

(defn test-reg
  [n1 n2]
  (println "\nfungp :: Functional Genetic Programming in Clojure")
  (println "Mike Vollmer, 2012")
  (println "Modified by Andrew Peklar and Eric Peck, Fall 2017")
  (println (str "Test inputs: " (vec in-list1)))
  (println (str "Test inputs: " (vec in-list2)))
  (println (str "Test outputs: " (vec out-list)))
  (println (str "Max generations: " (* n1 n2)))
  (println)
  (let [options {:iterations n1
                 :migrations n2
                 :num-islands 6
                 :population-size 40
                 :tournament-size 5
                 :mutation-probability 0.1
                 :max-depth 10
                 :terminals sample-parameters
                 :numbers number-literals
                 :fitness sample-fitness
                 :functions sample-functions
                 :report sample-report }
        [tree score] (rest (run-genetic-programming options))]
    (do (println "Done!")
        (sample-report tree score))))

