
(deftemplate SS
	( (Neg (70 1) (80 0))
	  (Zer (70 0) (80 1) (90 0))
	  (Pos (80 0) (90 1))
	)
)


(deftemplate PC
	( (Neg (10 1) (50 0))
	  (Zer (10 0) (50 1) (90 0))
	  (Pos (50 0) (90 1))
	)
)



(deftemplate FS
	( (BN (70 1) (75 0))
	  (N  (70 0) (75 1) (80 0))
	  (Z  (75 0) (80 1) (85 0))
	  (P  (80 0) (85 1) (90 0))
	  (BP (85 0) (90 1))
	)
)

			(defrule Get_Observations
				?i <- (initial-fact)
				 =>
				 (printout t "Enter angle: ")
				 (bind ?response (read))
				 (assert (crispSS ?response))
				 (printout t "Enter angle delta: ")
				 (bind ?response (read))
				 (assert (crispPC ?response))
				 (retract ?i))

(defrule fuzzify
	(crispSS ?a)
	(crispPC ?d)
	=>
	;(assert (SS (?a 0) (?a 1) (?a 0)))
	;(assert (PC (?d 0) (?d 1) (?d 0)))
	)

(defrule defuzzify1
(declare (salience -1))
?f <- (FS ?)
=>
(bind ?t (moment-defuzzify ?f))
(printout t "action--> " ?t crlf))







(defrule PP
	(SS Pos)
	(PC Pos)
	=>
	(assert (FS BP)))

	(defrule PZ
		(SS Pos)
		(PC Zer)
		=>   (assert (FS P)))

		(defrule PN
			(SS Pos)
			(PC Neg)
			=>
			(assert (FS Z)))


			(defrule ZP
				(SS Zer)
				(PC Pos)
				=>
				(assert (FS P)))

				(defrule ZZ
					(SS Zer)
					(PC Zer)
					=>
					(assert (FS Z)))

					(defrule ZN
						(SS Zer)
						(PC Neg)
						=>   (assert (FS N)))

						 (defrule NP
							 (SS Neg)
							 (PC Pos)
							 =>
							 (assert (FS Z)))

							  (defrule NZ
									 (SS Neg)
									 (PC Zer)
									 =>
									 (assert (FS N)))

									 (defrule NN
										 (SS Neg)
										 (PC Neg)
										 =>
										 (assert (FS BN)))
