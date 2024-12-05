# Open questions

- When merging the files and informations for computing the CProfiles, should we use the true or the weighted information?
- What is the range of momentum used for the tuning? (depends on particle type)
- How is the hybrid configuration made in WCSim?
    - Commented out line 386 of [WCSimEventAction.cc](http://WCSimEventAction.cc) (too much verbosity)
    - What happens when I do “hybrid=false”?
    - What happens when I set photo coverage for mPMT to 0?
- What range were used for fitting the charge PDF? and the polynomial degrees?
    - For HK,