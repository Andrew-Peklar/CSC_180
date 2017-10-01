;-------------------------------------------------
;Authors:	Andrew Peklar and Natalia Barakat
;Course :	CSC 180
;Section:	Natalia sec1, Andrew sec2\

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

;Facts:   Database lookup table for patient insurance accounts
;Fee  :   Contract/None
;Wait :   Preventative-Basic-Major
;Paid :   amount insurance has paid in past year
;Share:   ShareB = Premium Contribution (%) for Basic coverage.
;	    :   ShareM = .............................Major coverage.
;Max  :   Maximum amount insurance is willing to pay now (instantaneous payment)
;DOB  :   date of birth stored as integer (similur to how SQL systems stored them back in '08)
;	    :   20170929 (today) - yyyymmdd (their birthday)
;Proph:   Teeth cleaning + flouride preventative (visits in a year)
;           ----------    -----------   ----------   ----------   -----   ------    ----     ----    ------   ------   ------   -----      --------     ----------   ----
;             Acc ID         SSN           DOB        Eff Date    Proph   Deduct    Max      Paid    ShareB   ShareM   WaitB    WaitM       Type        Insureance   Fee
;           ----------    -----------   ----------   ----------   -----   ------    ----     ----    ------   ------   ------   -----      --------     ----------   ----
(deffacts		Accounts
  (Insuree  123456789A    123-25-3456   20050315     06/01/2016     1       50      1250      100      20       50       6       12        Private         Cross     Cont  )
  (Insuree  IXN453L339    123-24-3456   20010325     06/01/2016     1       50      1250      100      20       50       6       12        Private         Cross     Cont  )
  (Insuree  ANX098H115    339-00-4567   19850101     03/10/2010     0       55      2000      460      10       30       0        0        Employer       Shield     None  )
  (Insuree  XGZ123A457    560-12-3456   19680406     01/01/2014     2       45      0        2000      15       45       6        9        Private         Cross     Cont  )
  (Insuree  PLZ870K090    909-33-7651   19700910     06/01/2016     0       50      1250      278      20       50       6       12        Private         Cross     Cont  )
  (Insuree  ERS293Q814    451-19-5771   19891202     05/15/2013     1       55      1500      720      25       55       0        0        Employer       Shield     None  )
);          ----------    -----------   ----------   ----------   -----   ------    ----     ----    ------   ------   -----    -----      --------     ---------    ----
;  Label        ?A            ?B            ?C           ?D         ?E      ?F       ?G       ?H       ?I       ?J      ?K        ?L          ?M            ?N        ?O

;Facts:   This database contains insurance information about dental services
;CF   :   Contract Fee rates for services
;NCF  :   Non-Contract Fee rates for services
;SSS  :   Standard Share Shield percentage (%)
;SSC  :   ...............Cross................
;Cat  :   Insurance categories for patient services
;           ---   ---------------       --      ----    ---     ---    -----
;           Num   Service               CF       NCF    SSS     SSC    Cat.
;           ---   ---------------       --      ----    ---     ---    -----
(deffacts    Dental-Services
	(services  1    X-Ray                 30        30    100     100    Prev.  )
	(services  2    Examination           40        40    100     100    Prev.  )
	(services  3    Alloy-Filing          75       100    75       70    Basic  )
	(services  4    Install-Crown         700      800    50       45    Major  )
	(services  5    Composite-Fill        100      125    75       70    Basic  )
	(services  6    Routine-Extract       90       150    85       80    Basic  )
	(services  7    Fluoride-Cleaning     70        80    100     100    Prev.  )
	(services  8    Denture-Evaluation    1800    2000    70       65    Major  )
);           -    ------------------    ----    ----    ---    ----    -----
;           ?T      ?U                   ?V      ?W      ?X     ?Y      ?Z
;Facts: used to set input as facts
(deffacts patient-account-input
	(get-next-input)
);------------------------------------------------

(deffacts todays-date
	(date 20170928)
);------------------------------------------------

;Prompt: Start and ask for Insurance ID or SSN if ID not available
(defrule patient-ID-prompt
	?info <- (get-next-input)
	(date ?TD)
	=>
	(printout t " " crlf)
	(printout t "********************************************************************************************************* " crlf)
	(printout t "Disclaimer: Eligibility does not guarantee that the payment requested for the procedure will be provided. " crlf)
	(printout t "********************************************************************************************************* " crlf)
	(printout t " " crlf)
	(printout t "Hello! Welcome to the Dental Insurance Expert System " crlf)
	(printout t "Please have your information ready so that we may find your account" crlf)
	(printout t " " crlf)
	(printout t "Enter Insurance ID (use SSN if unknown): " )
		(bind ?A (read))
	(printout t "Enter Date of birth (yyyymmdd): " )
		(bind ?C (read))
	(printout t " " crlf)
	(assert (ID_or_SSN	?A))
	(assert (DOB-check	?C))
);-----------------------------------------------------------------------

;Prompt: If valid insurance ID is known
(defrule patient-has-ID
	?info <- 	(ID_or_SSN	?A)
	(DOB-check	?C)
	(Insuree	?A	?B	?C	?D	?E	?F	?G	?H	?I	?J	?K	?L	?M	?N ?O)
	=>
 	(printout t "Hello registered user ID Number " ?A " " crlf)
	(printout t " " crlf)
	(assert(patient-is ?A	?B	?C	?D	?E	?F	?G	?H	?I	?J	?K	?L	?M	?N ?O))
	(retract ?info)
);-----------------------------------------------------------------------

;Prompt: If SSN is entered instead of ID
(defrule patient-has-SSN
		?info <- (ID_or_SSN	?B)
		(DOB-check	?C)
		(Insuree	?A	?B	?C	?D	?E	?F	?G	?H	?I	?J	?K	?L	?M	?N ?O)
		=>
		(printout t "Hello registered user ID Number " ?A " " crlf)
		(printout t " " crlf)
		(assert(patient-is ?A	?B	?C	?D	?E	?F	?G	?H	?I	?J	?K	?L	?M	?N ?O))
		(retract ?info)
);-----------------------------------------------------------------------

;Prompt: Patient defines reason for visit and warnigs/disclaimers are issued
(defrule patient-service
	(get-next-input)
	(patient-is ?A	?B	?C	?D	?E	?F	?G	?H	?I	?J	?K	?L	?M	?N ?O)
	=>
	(printout t "What is the purpose for your visit today? " crlf)
	(printout t "Please enter the corrosponding service number " crlf )
		(printout t "  (1)  X-Ray " crlf)
		(printout t "  (2)  Examination " crlf)
		(printout t "  (3)  Alloy filing " crlf)
		(printout t "  (4)  Install a crown " crlf)
		(printout t "  (5)  Composite filing " crlf)
		(printout t "  (6)  Routine extraction " crlf)
		(printout t "  (7)  Cleaning + Flouride treatment " crlf)
		(printout t "  (8)  Upper and/or Lower denture eval " crlf)
	(printout t " " crlf)
	(printout t "Service number: " )
		(bind ?Service (read))
		(assert (selected-service ?Service))
);-----------------------------------------------------------------------

;ask if this is first visit this year
(defrule patient-q-first-visit
	(patient-is ?A	?B	?C	?D	?E	?F	?G	?H	?I	?J	?K	?L	?M	?N ?O)
	=>
	(printout t "Is this your first visit (yes/no): ")
		(bind ?YN (read))
	(assert (first-visit ?YN))
);-----------------------------------------------------------------------

(defrule prophy-Deny-Exceed-Private
	?kk <- (patient-is ?A	?B	?C	?D	2 ?F	?G	?H	?I	?J	?K	?L	?M	?N ?O)
	(selected-service 7)
	=>
	(printout t "You are ineligable for this procedure!" crlf )
	(printout t "Reason: have exceeded number of yearly cleanings " crlf )
	(retract ?kk)
)

(defrule prophy-Deny-Exceed-Employer
	?kk <- (patient-is ?A	?B	?C	?D	1 ?F	?G	?H	?I	?J	?K	?L	Employer	?N ?O)
	(selected-service 7)
	=>
	(printout t "You are ineligable for this procedure!" crlf )
	(printout t "Reason: have exceeded number of yearly cleanings " crlf )
	(retract ?kk)
)

(defrule Service-Crisp-Major-Cross
	(patient-is ?A	?B	?C	?D	?E ?F	?G	?H	?I	?J	?K	?L	?M	Cross ?O)
	(selected-service ?T)
	(services ?T ?U ?V ?W ?X ?Y Major)
	=>
	(assert (crispSSMajor ?Y))
	(assert (crispPC ?J))
)

(defrule Service-Crisp-Major-Shield
	(patient-is ?A	?B	?C	?D	?E ?F	?G	?H	?I	?J	?K	?L	?M	Shield ?O)
	(selected-service ?T)
	(services ?T ?U ?V ?W ?X ?Y Major)
	=>
	(assert (crispSSMajor ?X))
	(assert (crispPC ?J))
)

(defrule Service-Crisp-Basic-Cross
	(patient-is ?A	?B	?C	?D	?E ?F	?G	?H	?I	?J	?K	?L	?M	Cross ?O)
	(selected-service ?T)
	(services ?T ?U ?V ?W ?X ?Y Basic)
	=>
	(assert (crispSSBasic ?Y))
	(assert (crispPC ?I))
)

(defrule Service-Crisp-Basic-Shield
	(patient-is ?A	?B	?C	?D	?E ?F	?G	?H	?I	?J	?K	?L	?M	Shield ?O)
	(selected-service ?T)
	(services ?T ?U ?V ?W ?X ?Y Basic)
	=>
	(assert (crispSSBasic ?X))
	(assert (crispPC ?I))
)

(defrule fuzzify-Major
	(crispSSMajor ?a)
	(crispPC ?d)
	=>
	(assert (SSBasic (?a 0) (?a 1) (?a 0)))
	(assert (PC (?d 0) (?d 1) (?d 0)))
)

(defrule fuzzify-Basic
	(crispSSBasic ?a)
	(crispPC ?d)
	=>
	(assert (SSBasic (?a 0) (?a 1) (?a 0)))
	(assert (PC (?d 0) (?d 1) (?d 0)))
)


(defrule defuzzify-FSMajor
	(declare (salience -1))
	?f <- (FSMajor ?)
	=>
	(bind ?t (moment-defuzzify ?f))
	(printout t "action--> " ?t crlf)
)

(defrule defuzzify-FSBasic
	(declare (salience -1))
	?f <- (FSBasic ?)
	=>
	(bind ?t (moment-defuzzify ?f))
	(printout t "action--> " ?t crlf)
)

(defrule MPP
	(SSMajor Pos)
	(PC Pos)
	=>
	(assert (FSMajor BP)))

	(defrule MPZ
		(SSMajor Pos)
		(PC Zer)
		=>   (assert (FSMajor P)))

		(defrule MPN
			(SSMajor Pos)
			(PC Neg)
			=>
			(assert (FSMajor Z)))


			(defrule MZP
				(SSMajor Zer)
				(PC Pos)
				=>
				(assert (FSMajor P)))

				(defrule MZZ
					(SSMajor Zer)
					(PC Zer)
					=>
					(assert (FSMajor Z)))

					(defrule MZN
						(SSMajor Zer)
						(PC Neg)
						=>   (assert (FSMajor N)))

						 (defrule MNP
							 (SSMajor Neg)
							 (PC Pos)
							 =>
							 (assert (FSMajor Z)))

							  (defrule MNZ
									 (SSMajor Neg)
									 (PC Zer)
									 =>
									 (assert (FSMajor N)))

									 (defrule MNN
										 (SSMajor Neg)
										 (PC Neg)
										 =>
										 (assert (FSMajor BN)))

(defrule BPP
	(SSBasic Pos)
	(PC Pos)
	=>
	(assert (FSBasic BP)))

	(defrule BPZ
		(SSBasic Pos)
		(PC Zer)
		=>   (assert (FSBasic P)))

		(defrule BPN
			(SSBasic Pos)
			(PC Neg)
			=>
			(assert (FSBasic Z)))


			(defrule BZP
				(SSBasic Zer)
				(PC Pos)
				=>
				(assert (FSBasic P)))

				(defrule BZZ
					(SSBasic Zer)
					(PC Zer)
					=>
					(assert (FSBasic Z)))

					(defrule BZN
						(SSBasic Zer)
						(PC Neg)
						=>   (assert (FSBasic N)))

						 (defrule BNP
							 (SSBasic Neg)
							 (PC Pos)
							 =>
							 (assert (FSBasic Z)))

							  (defrule BNZ
									 (SSBasic Neg)
									 (PC Zer)
									 =>
									 (assert (FSBasic N)))

									 (defrule BNN
										 (SSBasic Neg)
										 (PC Neg)
										 =>
										 (assert (FSBasic BN)))
