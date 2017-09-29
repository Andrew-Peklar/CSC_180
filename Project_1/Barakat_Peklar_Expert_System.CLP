;-------------------------------------------------
;Authors:	Andrew Peklar and Natalia Barakat
;Course :	CSC 180
;Section:	Natalia sec1, Andrew sec2
;-------------------------------------------------

;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++Start of Facts++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
;           ----------    -----------   ----------   ----------   -----   ------    ----     ----    ------   ------   ----------------   --------    ----------    ----
;             Acc ID         SSN           DOB        Eff Date    Proph   Deduct    Max      Paid    ShareB   ShareM   Waiting P-B-M      Type        Insureance    Fee
;           ----------    -----------   ----------   ----------   -----   ------    ----     ----    ------   ------   ----------------   --------    ----------    ----
(deffacts		Accounts
  (Insuree  123456789A    123-25-3456   2008 0325    06/01/2016     1       50      1250      100      20       50      6 12               Private         Cross     Cont  )
  (Insuree  IXN453L339    123-24-3456   2007 0315    06/01/2016     1       50      1250      100      20       50      6 12               Private         Cross     Cont  )
  (Insuree  ANX098H115    339-00-4567   01/01/1985   03/10/2010     0       55      2000      460      10       30      0 0                Employer       Shield     None  )
  (Insuree  XGZ123A457    560-12-3456   04/06/1968   01/01/2014     2       45      0        2000      15       45      6month-9month      Private         Cross     Cont  )
  (Insuree  PLZ870K090    909-33-7651   09/10/1970   06/01/2016     0       50      1250      278      20       50      6month-1year       Private         Cross     Cont  )
  (Insuree  ERS293Q814    451-19-5771   12/02/1989   05/15/2013     1       55      1500      720      25       55      0-0                Employer       Shield     None  )
);          ----------    -----------   ----------   ----------   -----   ------    ----     ----    ------   ------   ---------------    ---------    ---------    ----
;  Label        ?A            ?B            ?C           ?D         ?E      ?F       ?G       ?H       ?I       ?J            ?K              ?L          ?M         ?N

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

;Facts: used to set input as facts
(deffacts patient-account-input
	(get-next-input)
);------------------------------------------------
;++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++End Of Facts+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++Start of User Input+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;Prompt: Start and ask for Insurance ID or SSN if ID not available
(defrule patient-ID-prompt
	?info <- (get-next-input)
	=>
	(printout t "Hello! Welcome to the Dental Insurance Expert System " crlf)
	(printout t "Please have your information ready so that we may find your account" crlf)
	(printout t "Enter Insurance ID (use SSN if unknown): " )
		(bind ?A (read))
	(assert (ID_or_SSN	?A))
);-----------------------------------------------------------------------

;Prompt: If valid insurance ID is known
(defrule patient-has-ID
	?info <- 	(ID_or_SSN	?A)
	(Insuree	?A	?B	?C	?D	?E	?F	?G	?H	?I	?J	?K	?L	?M	?N)
	=>
	(printout t "Enter Date of Birth (MM/DD/YYYY): ")
		(bind ?C (read))
	(assert(patient-is ?A	?B	?C	?D	?E	?F	?G	?H	?I	?J	?K	?L	?M	?N))
	(retract ?info)
);-----------------------------------------------------------------------

;Prompt: If SSN is entered instead of ID
(defrule patient-has-SSN
		?info <- (ID_or_SSN	?B)
		(Insuree	?A	?B	?C	?D	?E	?F	?G	?H	?I	?J	?K	?L	?M	?N)
		=>
		(printout t "Enter Date of Birth (MM DD YYYY): ")															;!!!!!!!!!!!
			(bind ?MM ?DD ?YY (read))
		(assert(patient-is ?A	?B	?C	?D	?E	?F	?G	?H	?I	?J	?K	?L	?M	?N))    ; make relational schema instead
		(retract ?info)
);-----------------------------------------------------------------------

;Prompt: Ask if first visit this year
(defrule patient-q-first-visit
	?info <- (get-next-input)
	=>
	(printout t "Is this your first visit (yes/no): ")
		(bind ?YN (read))
	(assert (first-visit ?YN))
);-----------------------------------------------------------------------

;TODO
(defrule patient-minor-or-adult

);-----------------------------------------------------------------------

;Prompt: Patient defines reason for visit and warnigs/disclaimers are issued
(defrule patient-service
	(get-next-input)
	=>
	(printout t "What is the purpose for your visit today? " crlf)
		(printout t "  <1> X-Ray " crlf)
		(printout t "  <2> Examination " crlf)
		(printout t "  <3> Alloy filing " crlf)
		(printout t "  <4> Install a crown " crlf)
		(printout t "  <5> Composite filing " crlf)
		(printout t "  <6> Routine extraction " crlf)
		(printout t "  <7> Cleaning + Flouride treatment " crlf)
		(printout t "  <8> Upper and/or Lower denture eval " crlf)
	(printout t "Please enter the corrosponding service number: " )
			(bind ?Service (read))
			(assert (selected-service ?Service))
	(printout t " " crlf)
	(printout t "Warning:    There is no 100% guarantee of service, regardless of coverage " crlf)
	(printout t "Disclaimer: Eligibility does not guaranty that the payment requested for the procedure will be provided " crlf)
);-----------------------------------------------------------------------


;Prompt:
(defrule patient-first-visit
	?info <- (first-visit	yes)
	=>
	(retract first-visit yes)
);-----------------------------------------------------------------------
