
(deftemplate SSMajor
	20 100 MinMax
	( (Neg (40 1) (60 0))
	  (Zer (40 0) (60 1) (80 0))
	  (Pos (60 0) (80 1))
	)
)

(deftemplate SSBasic
	60 100 MinMax
	( (Neg (70 1) (80 0))
	  (Zer (70 0) (80 1) (90 0))
	  (Pos (80 0) (90 1))
	)
)


(deftemplate PC
	0 100 MinMax
	( (Neg (10 1) (50 0))
	  (Zer (10 0) (50 1) (90 0))
	  (Pos (50 0) (90 1))
	)
)



(deftemplate FSMajor
	20 100 MinMax
	( (BN (40 1) (50 0))
	  (N  (40 0) (50 1) (60 0))
	  (Z  (50 0) (60 1) (70 0))
	  (P  (60 0) (70 1) (80 0))
	  (BP (70 0) (80 1))
	)
)

(deftemplate FSBasic
	60 100 MinMax
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
	(retract ?i)
)

(defrule fuzzify
	(crispSS ?a)
	(crispPC ?d)
	=>
	(assert (SSBasic (?a 0) (?a 1) (?a 0)))
	(assert (PC (?d 0) (?d 1) (?d 0))))

(defrule defuzzify1
(declare (salience -1))
?f <- (FSBasic ?)
=>
(bind ?t (moment-defuzzify ?f))
(printout t "action--> " ?t crlf))


(defrule PP
	(SSMajor Pos)
	(PC Pos)
	=>
	(assert (FSMajor BP)))

	(defrule PZ
		(SSMajor Pos)
		(PC Zer)
		=>   (assert (FSMajor P)))

		(defrule PN
			(SSMajor Pos)
			(PC Neg)
			=>
			(assert (FSMajor Z)))


			(defrule ZP
				(SSMajor Zer)
				(PC Pos)
				=>
				(assert (FSMajor P)))

				(defrule ZZ
					(SSMajor Zer)
					(PC Zer)
					=>
					(assert (FSMajor Z)))

					(defrule ZN
						(SSMajor Zer)
						(PC Neg)
						=>   (assert (FSMajor N)))

						 (defrule NP
							 (SSMajor Neg)
							 (PC Pos)
							 =>
							 (assert (FSMajor Z)))

							  (defrule NZ
									 (SSMajor Neg)
									 (PC Zer)
									 =>
									 (assert (FSMajor N)))

									 (defrule NN
										 (SSMajor Neg)
										 (PC Neg)
										 =>
										 (assert (FSMajor BN)))

(defrule PP
	(SSBasic Pos)
	(PC Pos)
	=>
	(assert (FSBasic BP)))

	(defrule PZ
		(SSBasic Pos)
		(PC Zer)
		=>   (assert (FSBasic P)))

		(defrule PN
			(SSBasic Pos)
			(PC Neg)
			=>
			(assert (FSBasic Z)))


			(defrule ZP
				(SSBasic Zer)
				(PC Pos)
				=>
				(assert (FSBasic P)))

				(defrule ZZ
					(SSBasic Zer)
					(PC Zer)
					=>
					(assert (FSBasic Z)))

					(defrule ZN
						(SSBasic Zer)
						(PC Neg)
						=>   (assert (FSBasic N)))

						 (defrule NP
							 (SSBasic Neg)
							 (PC Pos)
							 =>
							 (assert (FSBasic Z)))

							  (defrule NZ
									 (SSBasic Neg)
									 (PC Zer)
									 =>
									 (assert (FSBasic N)))

									 (defrule NN
										 (SSBasic Neg)
										 (PC Neg)
										 =>
										 (assert (FSBasic BN)))
