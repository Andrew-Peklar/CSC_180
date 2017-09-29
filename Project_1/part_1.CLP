;----------------------------------------------------
;Authors:	Andrew Peklar and Natalia Barakat
;Course :	CSC 180
;Section:	Natalia sec1, Andrew sec2
;----------------------------------------------------

;Facts:	connections between cities
(deffacts connections
	(from-to tahoe 		fairfield )
	(from-to fairfield	galt	  )
	(from-to chico		galt	  )
	(from-to galt		sacramento)
	(from-to sacramento	davis	  )
	(from-to chico		davis	  )
);----------------------------------------------------

;Facts: used to set input as facts
(deffacts user-inputs
	(get-next-input)
);----------------------------------------------------

;AskU: user for cities X and Y
(defrule get-inputs
	?f <- (get-next-input)
	=>
	(printout t "Enter start city:  ")
	(bind ?X (read))
	(printout t "Enter destination: ")
	(bind ?Y (read))
	(assert (travel ?X ?Y))
	(assert (city	?Y))
	(retract ?f)
);----------------------------------------------------

;Cond: when chico to !davis remove path
(defrule chico-to-not-davis
	?f <- (from-to chico davis)
	(travel chico ?y)
	(not (city davis))
	=>
	(retract ?f)
);----------------------------------------------------

;Goal: end when start 'x' is equal to end 'y'
(defrule goal-halt
	?f <- (travel ?y ?y)
	=>
	(retract ?f)
);----------------------------------------------------

;Rule: Do from-to, swap 'from' with 'to', print
(defrule swap-print
	?f <- (travel ?x ?y)
	(from-to ?x ?d)
	=>
	(printout t "drive from " ?x " to " ?d "  " crlf)
	(retract ?f)
	(assert (travel ?d ?y))
);----------------------------------------------------







	




