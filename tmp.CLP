;Facts: Database lookup table for patient insurance accounts
;             Acc ID          SSN          DOB        Eff Date    Proph   Deduct     Max     Paid    Share P/B/M    Waiting P-B-M     Type         Insureance Name        Fee
;           -----------   -----------   ----------   ----------   -----   ------    ----     ----    -----------    --------------    --------     -------------------    ----
(deffacts		Accounts
  (Insuree  1234567899A   123-25-3456   03/25/2007   06/01/2016     1       50      1250      100      0/20/50      0-6month-1year     Private     "Anthem Blue Cross"    Cont  )
  (Insuree  IXN453L3389   123-24-3456   03/15/2007   06/01/2016     1       50      1250      100      0/20/50      0-6month-1year     Private     "Anthem Blue Cross"    Cont  )
  (Insuree  ANX098H1155   339-00-4567   01/01/1985   03/10/2010     0       55      2000      460      0/10/30      0-0-0              Employer    "Blue Shield of CA"    None  )
  (Insuree  XGZ123A4567   560-12-3456   04/06/1968   01/01/2014     2       45         0     2000      0/15/45      0-6month-9month    Private     "Anthem Blue Cross"    Cont  )
  (Insuree  PLZ870K0990   909-33-7651   09/10/1970   06/01/2016     0       50      1250      278      0/20/50      0-6month-1year     Private     "Anthem Blue Cross"    Cont  )
  (Insuree  ERS293Q8134   451-19-5771   12/02/1989   05/15/2013     1       55      1500      720      0/25/55      0-0-0              Employer    "Blue Shield of CA"    None  )
);          -----------   -----------   ----------   ----------   -----   ------    ----     ----    -----------    ---------------    --------     ------------------    ----
;  Label        ?A            ?B            ?C           ?D         ?E       ?F      ?G       ?H         ?I               ?J             ?K                ?L              ?M
