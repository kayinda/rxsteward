#!/usr/bin/env python3
"""
Rx Steward — Antimicrobial Resistance & Stewardship Decision-Support Tool
======================================================================
A clinical decision-support tool for antimicrobial stewardship in
resource-limited settings.

Susceptibility data follows a 4-level geographic hierarchy with fallback:
  Uganda → East Africa → Africa → Global
Each organism displays the most granular level available.

Data sources:
  AMR — WHO GLASS 2025, WHO AWaRe 2023, CLSI/EUCAST, IDSA/ATS guidelines,
         Uganda tertiary hospital surveillance 2020-2023 (Muwanguzi et al. 2024),
         East African meta-analyses (Tadesse et al. 2017; PMC7409632; PMC12849242).
         Uganda Clinical Guidelines 2023 (MoH Uganda) — empiric therapy alignment.

Run:  pip install streamlit pandas plotly
      streamlit run rxsteward.py
"""

import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from collections import defaultdict

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                         1. DATA LAYER                                  ║
# ╚══════════════════════════════════════════════════════════════════════════╝

# ── 1A. ANTIBIOTIC DATABASE ──────────────────────────────────────────────
# AWaRe: A=Access, W=Watch, R=Reserve  |  Route: PO/IV/IM
ANTIBIOTICS = {
    # --- PENICILLINS ---
    "Amoxicillin":          {"class":"Penicillin","aware":"A","route":"PO","spectrum":"Gram+, some Gram−","notes":"First-line Access antibiotic for CAP, UTI, otitis media."},
    "Amoxicillin-clavulanate":{"class":"Penicillin+BLI","aware":"A","route":"PO/IV","spectrum":"Broad (incl. anaerobes, some ESBL−)","notes":"Beta-lactamase stable. Key Access drug."},
    "Ampicillin":           {"class":"Penicillin","aware":"A","route":"IV","spectrum":"Gram+, Listeria, Enterococcus","notes":"High resistance in E. coli (~80-100% East Africa)."},
    "Benzylpenicillin":     {"class":"Penicillin","aware":"A","route":"IV/IM","spectrum":"Streptococci, Treponema, Neisseria (declining)","notes":"Still first-line for GAS pharyngitis."},
    "Cloxacillin/Flucloxacillin":{"class":"Penicillin (anti-staph)","aware":"A","route":"PO/IV","spectrum":"MSSA","notes":"Not active vs MRSA."},
    "Piperacillin-tazobactam":{"class":"Penicillin+BLI","aware":"W","route":"IV","spectrum":"Broad incl. Pseudomonas, anaerobes","notes":"Hospital Watch antibiotic. Reserve for severe nosocomial."},
    # --- CEPHALOSPORINS ---
    "Cefalexin":            {"class":"Cephalosporin (1G)","aware":"A","route":"PO","spectrum":"Gram+, limited Gram−","notes":"Oral Access cephalosporin for SSTI, simple UTI."},
    "Cefuroxime":           {"class":"Cephalosporin (2G)","aware":"W","route":"PO/IV","spectrum":"Broader Gram−, H. influenzae","notes":"Watch. Surgical prophylaxis."},
    "Ceftriaxone":          {"class":"Cephalosporin (3G)","aware":"W","route":"IV/IM","spectrum":"Broad Gram−, meningeal penetration","notes":"Key Watch drug. Resistance rising 30-50% E. Africa (ESBL). First-line severe infections."},
    "Ceftazidime":          {"class":"Cephalosporin (3G anti-Pseudomonal)","aware":"W","route":"IV","spectrum":"Pseudomonas, Gram−","notes":"Anti-pseudomonal. Not active vs ESBL."},
    "Cefepime":             {"class":"Cephalosporin (4G)","aware":"W","route":"IV","spectrum":"Broad + Pseudomonas, AmpC stable","notes":"Watch. Stable to AmpC but not ESBL/carbapenemase."},
    "Ceftazidime-avibactam":{"class":"Cephalosporin+BLI","aware":"R","route":"IV","spectrum":"CRE (KPC, OXA-48), ESBL","notes":"RESERVE. Not active vs NDM (metallo-β-lactamase)."},
    # --- CARBAPENEMS ---
    "Meropenem":            {"class":"Carbapenem","aware":"W","route":"IV","spectrum":"Ultra-broad: ESBL, AmpC, anaerobes","notes":"Watch. Last reliable Gram− option before Reserve. Preserve."},
    "Imipenem-cilastatin":  {"class":"Carbapenem","aware":"W","route":"IV","spectrum":"Ultra-broad (not Stenotrophomonas)","notes":"Watch. Seizure risk at high dose."},
    "Ertapenem":            {"class":"Carbapenem","aware":"W","route":"IV/IM","spectrum":"Broad (not Pseudomonas, not Acinetobacter)","notes":"Once daily. Useful OPAT option."},
    # --- AMINOGLYCOSIDES ---
    "Gentamicin":           {"class":"Aminoglycoside","aware":"A","route":"IV/IM","spectrum":"Gram− aerobic, Staph synergy","notes":"Access. Nephro/ototoxic. Resistance 30-50% E. Africa Gram−."},
    "Amikacin":             {"class":"Aminoglycoside","aware":"A","route":"IV/IM","spectrum":"Broad Gram−, less resistance than gentamicin","notes":"Retains activity vs many gentamicin-R strains."},
    # --- FLUOROQUINOLONES ---
    "Ciprofloxacin":        {"class":"Fluoroquinolone","aware":"W","route":"PO/IV","spectrum":"Gram− incl. Pseudomonas, atypicals","notes":"Watch. Increasing resistance. Tendon/CNS risk."},
    "Levofloxacin":         {"class":"Fluoroquinolone","aware":"W","route":"PO/IV","spectrum":"Respiratory pathogens + Gram−","notes":"Watch. Respiratory fluoroquinolone."},
    "Moxifloxacin":         {"class":"Fluoroquinolone","aware":"W","route":"PO/IV","spectrum":"Broad respiratory + anaerobes","notes":"Watch. QTc prolongation risk."},
    # --- MACROLIDES ---
    "Azithromycin":         {"class":"Macrolide","aware":"W","route":"PO/IV","spectrum":"Atypicals, respiratory, STI (gonorrhoea declining)","notes":"Watch. Long tissue half-life."},
    "Erythromycin":         {"class":"Macrolide","aware":"A","route":"PO/IV","spectrum":"Gram+, atypicals","notes":"Access. GI side effects; prokinetic use."},
    # --- GLYCOPEPTIDES ---
    "Vancomycin":           {"class":"Glycopeptide","aware":"W","route":"IV (PO for C. diff)","spectrum":"MRSA, Enterococcus (not VanA-R)","notes":"Watch. TDM required (trough 15-20). Nephrotoxic."},
    "Teicoplanin":          {"class":"Glycopeptide","aware":"W","route":"IV/IM","spectrum":"MRSA (not VanA if vanA)","notes":"Watch. Less nephrotoxic than vancomycin; loading needed."},
    # --- OXAZOLIDINONES ---
    "Linezolid":            {"class":"Oxazolidinone","aware":"R","route":"PO/IV","spectrum":"MRSA, VRE, MDR-TB","notes":"RESERVE. Myelosuppression (monitor FBC weekly). 100% PO bioavailability."},
    # --- TETRACYCLINES ---
    "Doxycycline":          {"class":"Tetracycline","aware":"A","route":"PO","spectrum":"Atypicals, Rickettsia, Brucella, malaria (prophylaxis)","notes":"Access. Photosensitivity. Not <8 yr."},
    "Tigecycline":          {"class":"Glycylcycline","aware":"R","route":"IV","spectrum":"Broad incl. MRSA, ESBL, Acinetobacter (not Pseudomonas)","notes":"RESERVE. FDA black box: increased mortality."},
    # --- NITROIMIDAZOLES ---
    "Metronidazole":        {"class":"Nitroimidazole","aware":"A","route":"PO/IV","spectrum":"Anaerobes, C. difficile, protozoa","notes":"Access. First-line anaerobic cover. Avoid alcohol."},
    # --- SULFONAMIDES ---
    "Co-trimoxazole":       {"class":"Sulfonamide+Trimethoprim","aware":"A","route":"PO/IV","spectrum":"UTI, PJP, Nocardia, Toxoplasma","notes":"Access. High resistance in E. coli E. Africa (~70-80%). OI prophylaxis in HIV."},
    # --- POLYMYXINS ---
    "Colistin":             {"class":"Polymyxin","aware":"R","route":"IV/Inhaled","spectrum":"CRE, MDR Pseudomonas, MDR Acinetobacter","notes":"RESERVE last-resort. Nephrotoxic. mcr-1 plasmid-mediated R emerging."},
    # --- NITROFURANS ---
    "Nitrofurantoin":       {"class":"Nitrofuran","aware":"A","route":"PO","spectrum":"E. coli, Enterococcus (UTI only)","notes":"Access. Lower UTI only (no systemic levels). Low resistance."},
    # --- FOSFOMYCIN ---
    "Fosfomycin":           {"class":"Phosphonic acid","aware":"A","route":"PO (IV=R)","spectrum":"E. coli, Enterococcus (UTI)","notes":"Access (PO for UTI). RESERVE IV for systemic ESBL/CRE."},
    # --- ANTI-TB (selected) ---
    "Isoniazid":            {"class":"Anti-TB","aware":"A","route":"PO","spectrum":"M. tuberculosis","notes":"First-line TB. katG/inhA mutations → resistance."},
    "Rifampicin":           {"class":"Anti-TB","aware":"A","route":"PO/IV","spectrum":"M. tuberculosis, MRSA adjunct, prosthetic infections","notes":"First-line TB. rpoB mutations → MDR-TB. CYP inducer."},
    "Bedaquiline":          {"class":"Diarylquinoline","aware":"R","route":"PO","spectrum":"MDR/XDR-TB","notes":"RESERVE. ATP synthase inhibitor. QTc monitoring."},
    # --- OTHER ---
    "Clindamycin":          {"class":"Lincosamide","aware":"A","route":"PO/IV","spectrum":"Gram+, anaerobes, toxin suppression (GAS/S. aureus)","notes":"Access. Crosses bone. C. diff risk. D-test for iMLSB."},
    "Chloramphenicol":      {"class":"Amphenicol","aware":"A","route":"PO/IV","spectrum":"Broad, meningeal penetration","notes":"Access. Aplastic anaemia (rare but fatal). Still used E. Africa."},
    "Cefixime":             {"class":"Cephalosporin (3G oral)","aware":"W","route":"PO","spectrum":"Gram− (E. coli, Klebsiella, Gonorrhoea, Salmonella)","notes":"Watch. Oral 3G cephalosporin; widely used in Uganda (UCG). Alternative for gonorrhoea (400mg stat) and UTI. Not anti-pseudomonal. Susceptible to ESBL."},
    "Cefpodoxime":          {"class":"Cephalosporin (3G oral)","aware":"W","route":"PO","spectrum":"Gram−, respiratory pathogens","notes":"Watch. Oral 3G cephalosporin; step-down from IV ceftriaxone."},
    "Benzathine penicillin G":{"class":"Penicillin (long-acting)","aware":"A","route":"IM","spectrum":"Treponema pallidum, GAS (rheumatic fever prophylaxis)","notes":"Access. Single IM dose for primary/secondary syphilis. Monthly IM for rheumatic fever prophylaxis. No resistance in T. pallidum reported."},
    "Procaine penicillin":  {"class":"Penicillin (intermediate-acting)","aware":"A","route":"IM","spectrum":"Streptococci, syphilis, anthrax","notes":"Access. Used in combination with benzathine penicillin for neurosyphilis workup."},
    "Spectinomycin":        {"class":"Aminocyclitol","aware":"A","route":"IM","spectrum":"Neisseria gonorrhoeae","notes":"Access. IM alternative for gonorrhoea when ceftriaxone unavailable. Not effective for pharyngeal gonorrhoea."},
    "Fluconazole":          {"class":"Azole antifungal","aware":"A","route":"PO/IV","spectrum":"Candida spp., Cryptococcus neoformans","notes":"Not an antibiotic — included for HIV-associated infections. Maintenance/consolidation for cryptococcal meningitis. Oropharyngeal/vaginal candidiasis."},
    "Cefotaxime":           {"class":"Cephalosporin (3G)","aware":"W","route":"IV","spectrum":"Broad Gram−, meningeal penetration","notes":"Watch. Alternative to ceftriaxone; preferred in neonates (no calcium interaction). ESBL-susceptible."},
    "Cefiderocol":          {"class":"Siderophore cephalosporin","aware":"R","route":"IV","spectrum":"CRE incl NDM, MDR Pseudomonas, Acinetobacter","notes":"RESERVE. Active vs metallo-β-lactamases (NDM)."},
}

# ── 1B. BACTERIA DATABASE ────────────────────────────────────────────────
# susceptibility: % susceptible globally / East Africa where data available
BACTERIA = {
    "Escherichia coli": {
        "gram": "Gram-negative rod", "family": "Enterobacteriaceae",
        "diseases": ["UTI (most common cause)","Bacteraemia/sepsis","Intra-abdominal infection","Neonatal meningitis","Traveller\'s diarrhoea (ETEC)","HAP/VAP (rare)"],
        "resistance_mechanisms": ["ESBL (CTX-M-15 dominant E. Africa)","AmpC","Carbapenemase (NDM-1 emerging)","Fluoroquinolone target mutations"],
        "susceptibility": {
            "Global":       {"Ampicillin":35,"Amoxicillin-clavulanate":75,"Ceftriaxone":60,"Ceftazidime":60,"Meropenem":98,"Gentamicin":80,"Amikacin":95,"Ciprofloxacin":70,"Co-trimoxazole":60,"Nitrofurantoin":95,"Fosfomycin":95},
            "Africa":       {"Ampicillin":15,"Amoxicillin-clavulanate":55,"Ceftriaxone":30,"Ceftazidime":30,"Meropenem":92,"Gentamicin":60,"Amikacin":88,"Ciprofloxacin":45,"Co-trimoxazole":30,"Nitrofurantoin":92,"Fosfomycin":90},
            "East Africa":  {"Ampicillin":10,"Amoxicillin-clavulanate":50,"Ceftriaxone":30,"Ceftazidime":28,"Meropenem":90,"Gentamicin":55,"Amikacin":85,"Ciprofloxacin":40,"Co-trimoxazole":25,"Nitrofurantoin":90,"Fosfomycin":88},
            "Uganda":       {"Ampicillin":8,"Ceftriaxone":28,"Ciprofloxacin":40,"Chloramphenicol":65,"Meropenem":88,"Gentamicin":50,"Nitrofurantoin":90},
        },
        "sources": {"Global":"GLASS2025; LSHTM commentary (>40% R to 3GC globally)","Africa":"GLASS2025 AFRO (>70% E.coli R to 3GC)","East Africa":"Tadesse2017; EAMR_GENETICS2020","Uganda":"Muwanguzi et al. 2024 (PMC11053536): 28% ceftriaxone S overall; 40% cipro S; MRSA study"},
        "key_notes": "ESBL prevalence 40-65% in East African hospitals (CTX-M-15 dominant). Uganda tertiary data (2020-2023): only 28% ceftriaxone susceptible. Nitrofurantoin/fosfomycin retain activity for uncomplicated UTI."
    },
    "Klebsiella pneumoniae": {
        "gram": "Gram-negative rod", "family": "Enterobacteriaceae",
        "diseases": ["Hospital-acquired pneumonia","Bacteraemia/sepsis","UTI (complicated)","Liver abscess (hypervirulent K1/K2)","Surgical site infection","Neonatal sepsis (major E. Africa)"],
        "resistance_mechanisms": ["ESBL (CTX-M, SHV, TEM)","Carbapenemase (KPC, NDM-1, OXA-48)","Porin loss (OmpK35/36)","Colistin resistance (mgrB, mcr)"],
        "susceptibility": {
            "Global":       {"Ampicillin":0,"Amoxicillin-clavulanate":70,"Ceftriaxone":45,"Meropenem":92,"Gentamicin":75,"Amikacin":90,"Ciprofloxacin":70,"Co-trimoxazole":55,"Colistin":95},
            "Africa":       {"Ampicillin":0,"Amoxicillin-clavulanate":35,"Ceftriaxone":25,"Meropenem":80,"Gentamicin":50,"Amikacin":80,"Ciprofloxacin":40,"Co-trimoxazole":20,"Colistin":90},
            "East Africa":  {"Ampicillin":0,"Amoxicillin-clavulanate":35,"Ceftriaxone":25,"Meropenem":78,"Gentamicin":45,"Amikacin":78,"Ciprofloxacin":35,"Co-trimoxazole":18,"Colistin":88},
            "Uganda":       {"Ceftriaxone":42,"Ciprofloxacin":30,"Meropenem":75},
        },
        "sources": {"Global":"GLASS2025 (>55% K.pn R to 3GC)","Africa":"GLASS2025 AFRO; imipenem R increasing 15%/yr","East Africa":"Tadesse2017; EAMR_GENETICS2020","Uganda":"Muwanguzi 2024 (PMC11053536): 42% Enterobacterales ceftriaxone S; cipro & meropenem declining 55% & 47%"},
        "key_notes": "Leading cause of neonatal sepsis in SSA. Uganda: meropenem susceptibility declining rapidly (47% decrease 2020-2023). Carbapenem resistance (NDM-1) is a regional emergency. WHO Critical priority pathogen."
    },
    "Staphylococcus aureus": {
        "gram": "Gram-positive coccus (clusters)", "family": "Staphylococcaceae",
        "diseases": ["Skin/soft tissue infection (SSTI)","Bacteraemia/endocarditis","Osteomyelitis","Pneumonia (post-influenza, HAP)","Surgical site infection","Prosthetic joint infection","Toxic shock syndrome"],
        "resistance_mechanisms": ["mecA/mecC \u2192 MRSA (PBP2a)","blaZ \u2192 penicillinase","Fluoroquinolone target mutations","VISA/VRSA (rare)","iMLSB (erm genes)"],
        "susceptibility": {
            "Global":       {"Benzylpenicillin":15,"Cloxacillin/Flucloxacillin":72,"Vancomycin":99,"Linezolid":99,"Doxycycline":90,"Co-trimoxazole":90,"Clindamycin":75,"Gentamicin":85,"Ciprofloxacin":65,"Rifampicin":95},
            "Africa":       {"Benzylpenicillin":10,"Cloxacillin/Flucloxacillin":58,"Vancomycin":98,"Linezolid":98,"Doxycycline":82,"Co-trimoxazole":46,"Clindamycin":60,"Gentamicin":65,"Ciprofloxacin":50,"Rifampicin":90},
            "East Africa":  {"Benzylpenicillin":8,"Cloxacillin/Flucloxacillin":50,"Vancomycin":97,"Linezolid":97,"Doxycycline":80,"Co-trimoxazole":45,"Clindamycin":58,"Gentamicin":60,"Ciprofloxacin":48,"Rifampicin":88},
            "Uganda":       {"Cloxacillin/Flucloxacillin":44,"Vancomycin":95,"Linezolid":96,"Co-trimoxazole":40,"Gentamicin":55},
        },
        "sources": {"Global":"GLASS2025","Africa":"MRSA meta-analysis 2025 (PMC12849242): MRSA 42.2% Africa","East Africa":"MRSA E. Africa 20-55% variable","Uganda":"Muwanguzi 2024: MRSA 56% (44% methicillin S); gentamicin S decreased 37%"},
        "key_notes": "MRSA 56% in Ugandan tertiary hospitals (Muwanguzi 2024). Africa-wide MRSA ~42% (2025 meta-analysis). Vancomycin/linezolid remain reliable. CA-MRSA often TMP-SMX and doxycycline susceptible."
    },
    "Streptococcus pneumoniae": {
        "gram": "Gram-positive diplococcus", "family": "Streptococcaceae",
        "diseases": ["Community-acquired pneumonia (CAP)","Bacterial meningitis","Otitis media","Sinusitis","Bacteraemia"],
        "resistance_mechanisms": ["PBP alterations (penicillin non-susceptibility)","Macrolide efflux (mefA) and target modification (ermB)","Fluoroquinolone resistance (rare)"],
        "susceptibility": {
            "Global":       {"Benzylpenicillin":75,"Amoxicillin":90,"Ceftriaxone":95,"Azithromycin":70,"Levofloxacin":99,"Vancomycin":100,"Chloramphenicol":85,"Co-trimoxazole":55},
            "Africa":       {"Benzylpenicillin":68,"Amoxicillin":82,"Ceftriaxone":90,"Azithromycin":62,"Levofloxacin":98,"Vancomycin":100,"Chloramphenicol":78,"Co-trimoxazole":38},
            "East Africa":  {"Benzylpenicillin":70,"Amoxicillin":85,"Ceftriaxone":92,"Azithromycin":65,"Levofloxacin":98,"Vancomycin":100,"Chloramphenicol":80,"Co-trimoxazole":40},
        },
        "sources": {"Global":"GLASS2025","Africa":"GLASS2025 AFRO; PCV impact data","East Africa":"Tadesse2017; ESTIMATE"},
        "key_notes": "PCV13/PCV15 reducing IPD. Penicillin non-susceptibility is dose-dependent for non-meningitis (high-dose amoxicillin effective). Meningitis requires ceftriaxone."
    },
    "Enterococcus faecalis": {
        "gram": "Gram-positive coccus (chains)", "family": "Enterococcaceae",
        "diseases": ["UTI (complicated)","Endocarditis","Intra-abdominal infection","Bacteraemia (line-related)"],
        "resistance_mechanisms": ["Intrinsic cephalosporin resistance","vanA/vanB \u2192 VRE","High-level aminoglycoside resistance (HLAR)"],
        "susceptibility": {
            "Global":       {"Ampicillin":95,"Vancomycin":95,"Linezolid":99,"Nitrofurantoin":90,"Fosfomycin":85,"Gentamicin":60,"Doxycycline":70},
            "Africa":       {"Ampicillin":88,"Vancomycin":82,"Linezolid":97,"Nitrofurantoin":86,"Fosfomycin":78,"Gentamicin":48,"Doxycycline":62},
            "East Africa":  {"Ampicillin":85,"Vancomycin":78,"Linezolid":96,"Nitrofurantoin":85,"Fosfomycin":75,"Gentamicin":45,"Doxycycline":60},
            "Uganda":       {"Vancomycin":72,"Linezolid":95},
        },
        "sources": {"Global":"GLASS2025; MANDELL","Africa":"ESTIMATE from surveillance compilations","East Africa":"ESTIMATE","Uganda":"Muwanguzi 2024: Enterococcus vancomycin 72% S"},
        "key_notes": "Intrinsically resistant to all cephalosporins, clindamycin, TMP-SMX. Uganda: vancomycin susceptibility only 72% \u2014 VRE is a significant problem in tertiary hospitals."
    },
    "Enterococcus faecium": {
        "gram": "Gram-positive coccus (chains)", "family": "Enterococcaceae",
        "diseases": ["Nosocomial bacteraemia","UTI (complicated/catheter)","Endocarditis","Surgical site infection"],
        "resistance_mechanisms": ["Intrinsic ampicillin resistance","vanA/vanB \u2192 VRE","Linezolid resistance (rare but increasing)"],
        "susceptibility": {
            "Global":       {"Ampicillin":20,"Vancomycin":70,"Linezolid":98,"Daptomycin":95,"Tigecycline":95},
            "Africa":       {"Ampicillin":12,"Vancomycin":60,"Linezolid":95,"Daptomycin":92,"Tigecycline":92},
            "East Africa":  {"Ampicillin":10,"Vancomycin":55,"Linezolid":94,"Daptomycin":90,"Tigecycline":90},
        },
        "sources": {"Global":"GLASS2025; MANDELL","Africa":"ESTIMATE","East Africa":"ESTIMATE"},
        "key_notes": "Intrinsically ampicillin-resistant. VRE (vanA) rates increasing; WHO High-priority pathogen. Linezolid or daptomycin for VRE."
    },
    "Pseudomonas aeruginosa": {
        "gram": "Gram-negative rod (non-fermenter)", "family": "Pseudomonadaceae",
        "diseases": ["Hospital-acquired/ventilator-associated pneumonia","Burn wound infection","Chronic lung infection (CF/bronchiectasis)","UTI (catheter-associated)","Bacteraemia","Otitis externa (malignant)","Keratitis"],
        "resistance_mechanisms": ["AmpC hyperproduction","Porin loss (OprD \u2014 carbapenem R)","Efflux pumps (MexAB-OprM)","Metallo-\u03b2-lactamase (VIM, IMP, NDM)","Fluoroquinolone target mutations"],
        "susceptibility": {
            "Global":       {"Piperacillin-tazobactam":80,"Ceftazidime":80,"Cefepime":82,"Meropenem":85,"Gentamicin":80,"Amikacin":90,"Ciprofloxacin":72,"Colistin":95,"Cefiderocol":95},
            "Africa":       {"Piperacillin-tazobactam":68,"Ceftazidime":68,"Cefepime":70,"Meropenem":75,"Gentamicin":62,"Amikacin":84,"Ciprofloxacin":55,"Colistin":92,"Cefiderocol":93},
            "East Africa":  {"Piperacillin-tazobactam":65,"Ceftazidime":65,"Cefepime":68,"Meropenem":72,"Gentamicin":60,"Amikacin":82,"Ciprofloxacin":55,"Colistin":90,"Cefiderocol":92},
        },
        "sources": {"Global":"GLASS2025","Africa":"GLASS2025 AFRO + ESTIMATE","East Africa":"Tadesse2017 + ESTIMATE"},
        "key_notes": "Intrinsically resistant to many antibiotics. WHO Critical priority. Anti-pseudomonal agents required: pip-tazo, ceftazidime, cefepime, meropenem, aminoglycosides, ciprofloxacin."
    },
    "Acinetobacter baumannii": {
        "gram": "Gram-negative coccobacillus (non-fermenter)", "family": "Moraxellaceae",
        "diseases": ["VAP/HAP (ICU)","Wound infection (trauma/burns)","Bacteraemia","UTI (catheter)","Meningitis (post-neurosurgical)"],
        "resistance_mechanisms": ["OXA carbapenemases (OXA-23 dominant E. Africa)","NDM/VIM metallo-\u03b2-lactamase","Efflux pumps","Porin modifications","Pan-drug resistance emerging"],
        "susceptibility": {
            "Global":       {"Meropenem":50,"Amikacin":60,"Colistin":90,"Tigecycline":70,"Cefiderocol":80,"Co-trimoxazole":50},
            "Africa":       {"Meropenem":35,"Amikacin":48,"Colistin":86,"Tigecycline":62,"Cefiderocol":76,"Co-trimoxazole":38},
            "East Africa":  {"Meropenem":30,"Amikacin":45,"Colistin":85,"Tigecycline":60,"Cefiderocol":75,"Co-trimoxazole":35},
        },
        "sources": {"Global":"GLASS2025","Africa":"GLASS2025 AFRO (carbapenem R increasing 15%/yr fastest in Africa)","East Africa":"EAMR_GENETICS2020 + ESTIMATE"},
        "key_notes": "WHO #1 Critical priority pathogen (carbapenem-resistant). Extremely resistant in E. African ICUs. Colistin often last resort. Sulbactam combinations under investigation."
    },
    "Salmonella typhi": {
        "gram": "Gram-negative rod", "family": "Enterobacteriaceae",
        "diseases": ["Typhoid fever","Bacteraemia","Osteomyelitis (sickle cell)"],
        "resistance_mechanisms": ["MDR (Amp-Chlor-TMP)","Fluoroquinolone DCS (gyrA mutations)","ESBL (XDR H58 lineage)","Azithromycin resistance (emerging)"],
        "susceptibility": {
            "Global":       {"Ampicillin":55,"Chloramphenicol":55,"Co-trimoxazole":55,"Ciprofloxacin":50,"Azithromycin":90,"Ceftriaxone":95,"Meropenem":100},
            "Africa":       {"Ampicillin":48,"Chloramphenicol":48,"Co-trimoxazole":48,"Ciprofloxacin":42,"Azithromycin":88,"Ceftriaxone":93,"Meropenem":100},
            "East Africa":  {"Ampicillin":45,"Chloramphenicol":45,"Co-trimoxazole":45,"Ciprofloxacin":40,"Azithromycin":88,"Ceftriaxone":92,"Meropenem":100},
        },
        "sources": {"Global":"GLASS2025","Africa":"GLASS2025 + Andrews et al. Lancet Global Health 2023","East Africa":"ESTIMATE + typhoid surveillance"},
        "key_notes": "XDR S. Typhi (H58/4.3.1 lineage) spreading. Ceftriaxone for severe; azithromycin for uncomplicated. TCV vaccination recommended."
    },
    "Non-typhoidal Salmonella (iNTS)": {
        "gram": "Gram-negative rod", "family": "Enterobacteriaceae",
        "diseases": ["Invasive bacteraemia (HIV/malaria/malnutrition)","Gastroenteritis","Meningitis (children, SSA)","Osteomyelitis"],
        "resistance_mechanisms": ["MDR ST313 lineage (adapted to bloodstream, SSA)","ESBL","Fluoroquinolone resistance"],
        "susceptibility": {
            "Global":       {"Ampicillin":60,"Ceftriaxone":90,"Ciprofloxacin":85,"Azithromycin":85,"Meropenem":99,"Co-trimoxazole":50},
            "Africa":       {"Ampicillin":35,"Ceftriaxone":82,"Ciprofloxacin":72,"Azithromycin":82,"Meropenem":98,"Co-trimoxazole":28},
            "East Africa":  {"Ampicillin":30,"Ceftriaxone":80,"Ciprofloxacin":70,"Azithromycin":80,"Meropenem":98,"Co-trimoxazole":25},
        },
        "sources": {"Global":"GLASS2025","Africa":"iNTS ST313 studies + GLASS AFRO","East Africa":"Tadesse2017 + ESTIMATE"},
        "key_notes": "iNTS ST313 is a major killer in SSA, especially in HIV/malaria/malnourished children. MDR common. Blood culture essential \u2014 clinical presentation mimics malaria."
    },
    "Neisseria gonorrhoeae": {
        "gram": "Gram-negative diplococcus", "family": "Neisseriaceae",
        "diseases": ["Urethritis","Cervicitis","PID","Epididymo-orchitis","Disseminated gonococcal infection","Ophthalmia neonatorum"],
        "resistance_mechanisms": ["PBP modifications (penA mosaic)","tetM acquisition","Fluoroquinolone resistance (gyrA)","Azithromycin resistance (23S rRNA mutations)","Ceftriaxone MIC creep"],
        "susceptibility": {
            "Global":       {"Ceftriaxone":95,"Azithromycin":80,"Ciprofloxacin":30,"Spectinomycin":95,"Gentamicin":90},
            "Africa":       {"Ceftriaxone":93,"Azithromycin":76,"Ciprofloxacin":22,"Spectinomycin":93,"Gentamicin":88},
            "East Africa":  {"Ceftriaxone":93,"Azithromycin":75,"Ciprofloxacin":20,"Spectinomycin":92,"Gentamicin":88},
        },
        "sources": {"Global":"GLASS2025 (gonorrhoea section)","Africa":"GLASS2025 AFRO","East Africa":"WHO GASP + ESTIMATE"},
        "key_notes": "Dual therapy (ceftriaxone 500mg IM + azithromycin 1g PO) standard but azithromycin resistance threatens this. WHO priority pathogen."
    },
    "Neisseria meningitidis": {
        "gram": "Gram-negative diplococcus", "family": "Neisseriaceae",
        "diseases": ["Bacterial meningitis","Meningococcaemia/sepsis","Waterhouse-Friderichsen syndrome"],
        "resistance_mechanisms": ["Penicillin reduced susceptibility (PBP2 alterations)","Chloramphenicol resistance (rare)"],
        "susceptibility": {
            "Global":       {"Benzylpenicillin":85,"Ceftriaxone":99,"Chloramphenicol":95,"Ciprofloxacin":98,"Meropenem":100,"Rifampicin":99},
            "Africa":       {"Benzylpenicillin":82,"Ceftriaxone":99,"Chloramphenicol":92,"Ciprofloxacin":97,"Meropenem":100,"Rifampicin":98},
        },
        "sources": {"Global":"GLASS2025; MANDELL","Africa":"WHO meningitis belt surveillance + ESTIMATE"},
        "key_notes": "Meningitis belt of Africa (serogroup A \u2014 MenAfriVac impact dramatic; W and X emerging). Ceftriaxone first-line for suspected meningitis."
    },
    "Mycobacterium tuberculosis": {
        "gram": "Acid-fast bacillus", "family": "Mycobacteriaceae",
        "diseases": ["Pulmonary TB","TB meningitis","Miliary TB","Skeletal TB","Genitourinary TB","Pericardial TB","Lymph node TB"],
        "resistance_mechanisms": ["katG S315T (INH high-level R)","inhA promoter (INH low-level R, ethionamide cross-R)","rpoB S531L (RIF R \u2014 most common)","embB M306V (EMB R)","pncA (PZA R)","gyrA/B (FQ R)","rrs (aminoglycoside R)"],
        "susceptibility": {
            "Global":       {"Isoniazid":87,"Rifampicin":95,"Ethambutol":95,"Pyrazinamide":96,"Levofloxacin":97,"Amikacin":98,"Bedaquiline":99,"Linezolid":99},
            "Africa":       {"Isoniazid":86,"Rifampicin":94,"Ethambutol":94,"Pyrazinamide":95,"Levofloxacin":97,"Amikacin":98,"Bedaquiline":99,"Linezolid":99},
            "East Africa":  {"Isoniazid":85,"Rifampicin":93,"Ethambutol":93,"Pyrazinamide":94,"Levofloxacin":96,"Amikacin":97,"Bedaquiline":99,"Linezolid":99},
            "Uganda":       {"Isoniazid":85,"Rifampicin":94},
        },
        "sources": {"Global":"WHO Global TB Report 2024","Africa":"WHO Africa TB Report","East Africa":"Tadesse2017 + NTLP data","Uganda":"Uganda NTLP Annual Report 2023; MDR-TB ~1.6% new cases"},
        "key_notes": "MDR-TB (INH+RIF R): ~3.3% new, ~18% retreatment globally. GeneXpert detects RIF-R in 2 hours. Uganda MDR-TB: ~1.6% new cases. Bedaquiline-based short regimens (BPaL) transforming treatment."
    },
    "Haemophilus influenzae": {
        "gram": "Gram-negative coccobacillus", "family": "Pasteurellaceae",
        "diseases": ["Otitis media","Sinusitis","CAP","Exacerbation of COPD","Meningitis (type b, declining post-Hib vaccine)","Epiglottitis"],
        "resistance_mechanisms": ["Beta-lactamase (TEM-1) \u2014 20-40%","BLNAR (PBP3 mutation)"],
        "susceptibility": {
            "Global":       {"Amoxicillin":75,"Amoxicillin-clavulanate":98,"Ceftriaxone":99,"Azithromycin":99,"Ciprofloxacin":99,"Chloramphenicol":95,"Co-trimoxazole":65},
            "Africa":       {"Amoxicillin":62,"Amoxicillin-clavulanate":94,"Ceftriaxone":98,"Azithromycin":96,"Ciprofloxacin":98,"Chloramphenicol":88,"Co-trimoxazole":48},
            "East Africa":  {"Amoxicillin":65,"Amoxicillin-clavulanate":95,"Ceftriaxone":98,"Azithromycin":97,"Ciprofloxacin":98,"Chloramphenicol":90,"Co-trimoxazole":50},
        },
        "sources": {"Global":"GLASS2025; MANDELL","Africa":"ESTIMATE","East Africa":"Tadesse2017 + ESTIMATE"},
        "key_notes": "Hib vaccine dramatically reduced invasive disease. Non-typeable H. influenzae now dominant. Beta-lactamase production 20-40% \u2014 amoxicillin-clavulanate covers."
    },
    "Group A Streptococcus (S. pyogenes)": {
        "gram": "Gram-positive coccus (chains)", "family": "Streptococcaceae",
        "diseases": ["Pharyngitis/tonsillitis","Cellulitis/erysipelas","Necrotising fasciitis","Impetigo","Rheumatic fever (sequela)","Post-streptococcal GN","Scarlet fever","Puerperal sepsis"],
        "resistance_mechanisms": ["Macrolide resistance (ermB, mefA) \u2014 5-30%","Tetracycline resistance \u2014 20-40%","No penicillin resistance reported"],
        "susceptibility": {
            "Global":       {"Benzylpenicillin":100,"Amoxicillin":100,"Ceftriaxone":100,"Azithromycin":80,"Clindamycin":90,"Vancomycin":100,"Doxycycline":70},
            "Africa":       {"Benzylpenicillin":100,"Amoxicillin":100,"Ceftriaxone":100,"Azithromycin":74,"Clindamycin":84,"Vancomycin":100,"Doxycycline":62},
            "East Africa":  {"Benzylpenicillin":100,"Amoxicillin":100,"Ceftriaxone":100,"Azithromycin":75,"Clindamycin":85,"Vancomycin":100,"Doxycycline":65},
        },
        "sources": {"Global":"MANDELL; universal penicillin susceptibility documented","Africa":"ESTIMATE","East Africa":"ESTIMATE"},
        "key_notes": "Penicillin remains 100% effective \u2014 the gold standard. Clindamycin added for toxin suppression in severe invasive GAS (necrotising fasciitis, toxic shock)."
    },
    "Group B Streptococcus (S. agalactiae)": {
        "gram": "Gram-positive coccus (chains)", "family": "Streptococcaceae",
        "diseases": ["Neonatal early-onset sepsis","Neonatal meningitis","Chorioamnionitis","UTI in pregnancy","Bacteraemia (elderly/diabetic)"],
        "resistance_mechanisms": ["Macrolide resistance (erm genes)","Clindamycin resistance","No penicillin resistance"],
        "susceptibility": {
            "Global":       {"Benzylpenicillin":100,"Ampicillin":100,"Ceftriaxone":100,"Vancomycin":100,"Clindamycin":80,"Azithromycin":70},
            "Africa":       {"Benzylpenicillin":100,"Ampicillin":100,"Ceftriaxone":100,"Vancomycin":100,"Clindamycin":74,"Azithromycin":64},
        },
        "sources": {"Global":"MANDELL; universal penicillin susceptibility","Africa":"ESTIMATE from African GBS studies"},
        "key_notes": "Universal penicillin susceptibility. IAP with penicillin G/ampicillin for GBS-colonised mothers. Clindamycin if penicillin-allergic AND D-test negative."
    },
    "Clostridioides difficile": {
        "gram": "Gram-positive anaerobic rod (spore-forming)", "family": "Peptostreptococcaceae",
        "diseases": ["Antibiotic-associated diarrhoea","Pseudomembranous colitis","Toxic megacolon","Recurrent CDI"],
        "resistance_mechanisms": ["Fluoroquinolone resistance (epidemic ribotype 027)","Metronidazole treatment failure (~20%)","Fidaxomicin resistance (rare)"],
        "susceptibility": {
            "Global":       {"Vancomycin":98,"Metronidazole":80,"Fidaxomicin":99},
        },
        "sources": {"Global":"IDSA/SHEA CDI guidelines 2021; MANDELL"},
        "key_notes": "Oral vancomycin (not IV) is first-line for CDI (IDSA 2021). Metronidazole only if vancomycin unavailable. Fidaxomicin reduces recurrence."
    },
    "Cryptococcus neoformans": {
        "gram": "Encapsulated yeast (not bacterial \u2014 included for HIV context)", "family": "Cryptococcaceae",
        "diseases": ["Cryptococcal meningitis (HIV, CD4<100)","Pulmonary cryptococcosis","Disseminated disease"],
        "resistance_mechanisms": ["Fluconazole heteroresistance","Amphotericin B tolerance"],
        "susceptibility": {
            "Global":       {"Amphotericin B":95,"Flucytosine":90,"Fluconazole":85},
            "Africa":       {"Amphotericin B":93,"Flucytosine":88,"Fluconazole":80},
            "East Africa":  {"Amphotericin B":93,"Flucytosine":88,"Fluconazole":80},
        },
        "sources": {"Global":"WHO Crypto guidelines 2022; Rajasingham et al. Lancet ID 2022","Africa":"ESTIMATE from African crypto studies","East Africa":"ESTIMATE"},
        "key_notes": "Leading cause of meningitis death in PLHIV in SSA. WHO 2022: 1 week AmB + flucytosine induction \u2192 fluconazole consolidation/maintenance. CrAg screening at CD4<200 saves lives."
    },
    "Treponema pallidum": {
        "gram": "Spirochaete (not Gram-stainable)", "family": "Spirochaetaceae",
        "diseases": ["Primary syphilis (chancre)","Secondary syphilis (rash, condylomata)","Tertiary syphilis (gumma, cardiovascular)","Neurosyphilis","Congenital syphilis"],
        "resistance_mechanisms": ["No penicillin resistance documented","Azithromycin resistance (23S rRNA A2058G \u2014 30-90%)","Tetracycline resistance (rare)"],
        "susceptibility": {"Global":{"Benzathine penicillin G":100,"Doxycycline":95,"Ceftriaxone":99,"Azithromycin":70},"East Africa":{"Benzathine penicillin G":100,"Doxycycline":95,"Ceftriaxone":99,"Azithromycin":60}},
        "sources": {"Global":"WHO STI guidelines 2021; CDC 2021","East Africa":"ESTIMATE"},
        "key_notes": "Benzathine penicillin G 100% effective \u2014 no resistance ever. Azithromycin NOT recommended (resistance). Screen all pregnant women."
    },
    "Vibrio cholerae": {
        "gram": "Gram-negative curved rod", "family": "Vibrionaceae",
        "diseases": ["Cholera (acute watery diarrhoea)","Dehydration/hypovolaemic shock"],
        "resistance_mechanisms": ["TMP-SMX resistance (widespread)","Tetracycline resistance (some strains)","Fluoroquinolone resistance (emerging)"],
        "susceptibility": {"Global":{"Doxycycline":90,"Azithromycin":95,"Ciprofloxacin":85,"Co-trimoxazole":40},"East Africa":{"Doxycycline":85,"Azithromycin":92,"Ciprofloxacin":80,"Co-trimoxazole":30}},
        "sources": {"Global":"WHO cholera guidelines 2023","East Africa":"ESTIMATE"},
        "key_notes": "ORS is the PRIMARY treatment. Antibiotics adjunctive. Doxycycline 300mg single dose or azithromycin 1g stat. Notifiable disease."
    },
    "Shigella spp.": {
        "gram": "Gram-negative rod", "family": "Enterobacteriaceae",
        "diseases": ["Bacillary dysentery (bloody diarrhoea)","Shigellosis","HUS (S. dysenteriae type 1)"],
        "resistance_mechanisms": ["Ampicillin/TMP-SMX resistance (widespread)","Fluoroquinolone resistance (emerging)"],
        "susceptibility": {"Global":{"Ciprofloxacin":75,"Azithromycin":85,"Ceftriaxone":95,"Ampicillin":30,"Co-trimoxazole":30},"East Africa":{"Ciprofloxacin":65,"Azithromycin":82,"Ceftriaxone":92,"Ampicillin":20,"Co-trimoxazole":20}},
        "sources": {"Global":"GLASS 2025","East Africa":"Tadesse2017 + ESTIMATE"},
        "key_notes": "Ciprofloxacin or azithromycin for confirmed shigellosis. Rising FQ resistance. Ceftriaxone IV for severe. Ampicillin/TMP-SMX no longer reliable."
    },
    "Helicobacter pylori": {
        "gram": "Gram-negative spiral rod", "family": "Helicobacteraceae",
        "diseases": ["Peptic ulcer disease","Gastritis","Gastric MALT lymphoma","Gastric adenocarcinoma (risk factor)"],
        "resistance_mechanisms": ["Clarithromycin resistance (15-40%)","Metronidazole resistance (40-90%)","Amoxicillin resistance (rare)"],
        "susceptibility": {"Global":{"Amoxicillin":97,"Metronidazole":40,"Clarithromycin":70,"Tetracycline":95},"Africa":{"Amoxicillin":95,"Metronidazole":25,"Clarithromycin":65,"Tetracycline":92}},
        "sources": {"Global":"Maastricht VI consensus 2022","Africa":"ESTIMATE"},
        "key_notes": "Triple therapy: PPI + amoxicillin + clarithromycin 14d. Metronidazole resistance very high in Africa (>60%). Consider bismuth quadruple if clarithromycin resistance >15%."
    },
    "Bordetella pertussis": {
        "gram": "Gram-negative coccobacillus", "family": "Alcaligenaceae",
        "diseases": ["Pertussis (whooping cough)","Neonatal pertussis"],
        "resistance_mechanisms": ["Macrolide resistance (very rare)"],
        "susceptibility": {"Global":{"Azithromycin":99,"Erythromycin":99,"Co-trimoxazole":90}},
        "sources": {"Global":"CDC pertussis guidelines; WHO"},
        "key_notes": "Macrolides are the ONLY effective treatment. DPT vaccination is primary prevention."
    },
    "Chlamydia trachomatis": {
        "gram": "Obligate intracellular bacterium", "family": "Chlamydiaceae",
        "diseases": ["Urethritis/cervicitis","PID","Neonatal conjunctivitis","Trachoma","LGV"],
        "resistance_mechanisms": ["No clinically significant resistance documented"],
        "susceptibility": {"Global":{"Doxycycline":99,"Azithromycin":95,"Levofloxacin":95}},
        "sources": {"Global":"WHO STI guidelines 2021; CDC 2021"},
        "key_notes": "Doxycycline 100mg BD \u00d7 7d now preferred over azithromycin 1g stat (higher cure rate). Screen pregnant women. Treat partners."
    },
    "Mycoplasma genitalium": {
        "gram": "Cell-wall-deficient bacterium", "family": "Mycoplasmataceae",
        "diseases": ["Non-gonococcal urethritis","Cervicitis","PID"],
        "resistance_mechanisms": ["Macrolide resistance (40-80%)","Intrinsic beta-lactam resistance (no cell wall)"],
        "susceptibility": {"Global":{"Doxycycline":40,"Azithromycin":50,"Moxifloxacin":90}},
        "sources": {"Global":"Australasian/European M. gen guidelines 2022"},
        "key_notes": "Resistance-guided therapy preferred. Doxycycline first to reduce load, then azithromycin or moxifloxacin. Intrinsically resistant to ALL beta-lactams."
    },
}

# ── 1C. CLINICAL SYNDROME → EMPIRIC THERAPY ─────────────────────────────
SYNDROMES = {
    "Uncomplicated UTI (cystitis)": {
        "setting": "Community",
        "likely_pathogens": ["Escherichia coli","Klebsiella pneumoniae","Staphylococcus saprophyticus","Enterococcus faecalis"],
        "first_line": [("Nitrofurantoin","100mg PO BD × 5d","A"),("Fosfomycin","3g PO single dose","A")],
        "second_line": [("Amoxicillin-clavulanate","625mg PO TDS × 5d","A"),("Cefalexin","500mg PO QDS × 5d","A")],
        "avoid": "Co-trimoxazole (>20% E. coli resistance in E. Africa); Ciprofloxacin (collateral damage — preserve for complicated UTI).",
        "notes": "Urine culture before antibiotics if recurrent/complicated. High co-trimoxazole resistance in E. Africa limits its use."
    },
    "Complicated UTI / Pyelonephritis": {
        "setting": "Community / Hospital",
        "likely_pathogens": ["Escherichia coli","Klebsiella pneumoniae","Proteus mirabilis","Pseudomonas aeruginosa (catheter)","Enterococcus spp."],
        "first_line": [("Ceftriaxone","1-2g IV OD","W"),("Gentamicin","5-7mg/kg IV OD","A")],
        "second_line": [("Ciprofloxacin","500mg PO BD (if susceptible)","W"),("Amoxicillin-clavulanate","1.2g IV TDS","A")],
        "avoid": "Nitrofurantoin (no systemic levels). Empiric ciprofloxacin if local FQ resistance >20%.",
        "notes": "Blood cultures + urine culture BEFORE antibiotics. De-escalate at 48-72h based on susceptibility."
    },
    "Community-acquired pneumonia (CAP)": {
        "setting": "Community / Hospital",
        "likely_pathogens": ["Streptococcus pneumoniae","Haemophilus influenzae","Mycoplasma pneumoniae","Legionella pneumophila","Staphylococcus aureus (post-influenza)"],
        "first_line": [("Amoxicillin","500mg-1g PO Q8h × 5d (mild/moderate, ambulatory)","A"),("Benzylpenicillin","2MU IV Q4-6h (severe, hospitalised) [UCG first-line]","A")],
        "second_line": [("Ceftriaxone","1g IV OD (if no response to penicillin at 48h) [UCG second-line]","W"),("Doxycycline","100mg PO BD × 7-10d (atypical/penicillin allergy) [UCG alternative]","A"),("Ceftriaxone + Azithromycin","IV (severe/ICU — ATS/IDSA guideline)","W+W")],
        "avoid": "Fluoroquinolone monotherapy as first-line (preserve for penicillin allergy/confirmed atypical). Co-trimoxazole (poor pneumococcal activity).",
        "notes": "CURB-65 score guides setting. UCG 2023: benzylpenicillin is first-line for severe CAP (Access antibiotic); ceftriaxone is second-line if no response at 48h. ATS/IDSA guidelines use ceftriaxone + azithromycin for severe — both are listed. Add PJP cover (co-trimoxazole) if HIV + CD4<200. Add TB workup if HIV or >2 weeks symptoms."
    },
    "Hospital-acquired / Ventilator-associated pneumonia (HAP/VAP)": {
        "setting": "Hospital / ICU",
        "likely_pathogens": ["Klebsiella pneumoniae","Pseudomonas aeruginosa","Acinetobacter baumannii","Staphylococcus aureus (MRSA)","Escherichia coli"],
        "first_line": [("Piperacillin-tazobactam","4.5g IV Q6h","W"),("Meropenem","1g IV Q8h (if ESBL/high resistance risk)","W")],
        "second_line": [("Cefepime + Amikacin","IV combination","W+A"),("Meropenem + Colistin","(if CRE/CRAB suspected)","W+R")],
        "avoid": "Narrow-spectrum agents alone. Monotherapy if MDR risk factors.",
        "notes": "Double Gram-negative cover if MDR risk factors (prior antibiotics, >5d hospitalisation, local ESBL/CRE >20%). Add vancomycin if MRSA risk. De-escalate aggressively."
    },
    "Bacterial meningitis (community-acquired)": {
        "setting": "Hospital (Emergency)",
        "likely_pathogens": ["Streptococcus pneumoniae","Neisseria meningitidis","Listeria monocytogenes (elderly/immunocompromised)","Group B Streptococcus (neonates)","E. coli (neonates)"],
        "first_line": [("Ceftriaxone","2g IV BD","W"),("+ Ampicillin","2g IV Q4h (if >50yr or immunocompromised — Listeria cover)","A")],
        "second_line": [("Meropenem","2g IV Q8h (penicillin allergy or resistant organisms)","W"),("Chloramphenicol","(if ceftriaxone unavailable — E. Africa resource setting)","A")],
        "avoid": "Any delay in antibiotics. Do NOT wait for CT/LP if clinical suspicion strong — treat immediately.",
        "notes": "Dexamethasone 0.15mg/kg Q6h × 4d BEFORE or WITH first dose of antibiotics (reduces mortality in pneumococcal meningitis). LP before antibiotics ONLY if no contraindication and no delay."
    },
    "Neonatal sepsis (early-onset)": {
        "setting": "Hospital (Neonatal unit)",
        "likely_pathogens": ["Group B Streptococcus","Escherichia coli","Listeria monocytogenes","Klebsiella pneumoniae (E. Africa)"],
        "first_line": [("Ampicillin + Gentamicin","IV","A+A")],
        "second_line": [("Ceftriaxone","(avoid in neonates with jaundice/calcium infusion)","W"),("Meropenem","(if ESBL suspected — E. Africa neonatal units)","W")],
        "avoid": "Ceftriaxone co-administered with calcium-containing IV solutions (risk of fatal precipitation).",
        "notes": "Klebsiella pneumoniae ESBL is the leading cause of neonatal sepsis deaths in many E. African units. Blood culture before antibiotics. Meropenem may be needed empirically in high ESBL-prevalence units."
    },
    "Skin and soft tissue infection (SSTI)": {
        "setting": "Community / Hospital",
        "likely_pathogens": ["Staphylococcus aureus","Group A Streptococcus","Anaerobes (diabetic foot, bite wounds)","Gram-negatives (diabetic foot)"],
        "first_line": [("Cloxacillin/Flucloxacillin","500mg PO QDS (mild, MSSA)","A"),("Amoxicillin-clavulanate","625mg PO TDS (bite/diabetic foot)","A")],
        "second_line": [("Clindamycin","450mg PO TDS (penicillin allergy/MRSA risk)","A"),("Vancomycin","IV (severe MRSA)","W")],
        "avoid": "Ciprofloxacin monotherapy for cellulitis (poor Gram+ activity).",
        "notes": "I&D for abscess is the primary treatment — antibiotics adjunctive. MRSA risk factors: prior MRSA, IVDU, incarceration, healthcare contact. Add anaerobic cover for diabetic foot/bite."
    },
    "Intra-abdominal infection": {
        "setting": "Hospital (Surgical)",
        "likely_pathogens": ["Escherichia coli","Klebsiella pneumoniae","Bacteroides fragilis","Enterococcus spp.","Pseudomonas aeruginosa (post-operative)"],
        "first_line": [("Ceftriaxone + Metronidazole","IV","W+A"),("Amoxicillin-clavulanate","1.2g IV TDS (mild-moderate)","A")],
        "second_line": [("Piperacillin-tazobactam","4.5g IV Q6h","W"),("Meropenem","(if ESBL risk)","W")],
        "avoid": "Source control failure — antibiotics alone won't treat undrained collections.",
        "notes": "Source control (surgery/drainage) is paramount. Duration: 4 days post-adequate source control (STOP-IT trial). Longer only if ongoing source."
    },
    "Infective endocarditis": {
        "setting": "Hospital",
        "likely_pathogens": ["Staphylococcus aureus","Viridans streptococci","Enterococcus faecalis","HACEK organisms","Coagulase-negative staphylococci (prosthetic valve)"],
        "first_line": [("Flucloxacillin + Gentamicin","IV (native valve, empiric)","A+A"),("Vancomycin + Gentamicin","(prosthetic valve or MRSA risk)","W+A")],
        "second_line": [("Ceftriaxone","2g IV OD × 4w (viridans strep, penicillin allergy)","W"),("Vancomycin + Rifampicin + Gentamicin","(prosthetic valve staph)","W+A+A")],
        "avoid": "Empiric treatment without blood cultures (3 sets from different sites before antibiotics).",
        "notes": "Modified Duke criteria for diagnosis. 4-6 weeks IV therapy. Early surgery consultation for heart failure, uncontrolled infection, large vegetations >10mm."
    },
    "Sepsis / Septic shock (undifferentiated)": {
        "setting": "Hospital / ICU",
        "likely_pathogens": ["Escherichia coli","Klebsiella pneumoniae","Staphylococcus aureus","Streptococcus pneumoniae","Pseudomonas aeruginosa"],
        "first_line": [("Gentamicin 7mg/kg IV OD + Cloxacillin 2g IV Q4-6h","(unknown focus, community-onset) [UCG first-line]","A+A"),("Ceftriaxone","2g IV OD (community-onset, alternative)","W")],
        "second_line": [("Meropenem","1g IV Q8h (hospital-onset / ESBL risk)","W"),("Piperacillin-tazobactam + Vancomycin","(broad cover + MRSA)","W+W"),("Meropenem + Amikacin + Vancomycin","(ICU, MDR risk)","W+A+W")],
        "avoid": "Any delay. Surviving Sepsis: antibiotics within 1 hour of recognition.",
        "notes": "Source identification critical. UCG 2023: gentamicin + cloxacillin is first-line for unknown-focus sepsis (both Access antibiotics). Ceftriaxone is the alternative. Blood cultures × 2 BEFORE antibiotics but do NOT delay. Lactate, fluid resuscitation, vasopressors per Surviving Sepsis 2021."
    },
    "Sexually transmitted infection (urethral/cervical discharge)": {
        "setting": "Community",
        "likely_pathogens": ["Neisseria gonorrhoeae","Chlamydia trachomatis","Mycoplasma genitalium","Trichomonas vaginalis"],
        "first_line": [("Ceftriaxone 500mg IM + Azithromycin 1g PO","single dose (syndromic)","W+W"),("Doxycycline","100mg PO BD × 7d (chlamydia)","A")],
        "second_line": [("Cefixime 400mg PO + Azithromycin","(if IM not available)","W+W")],
        "avoid": "Ciprofloxacin for gonorrhoea (>50% resistance in E. Africa). Azithromycin monotherapy (resistance rising).",
        "notes": "Syndromic management (WHO) when lab unavailable: treat for both gonorrhoea AND chlamydia simultaneously. Partner notification essential."
    },
    "Typhoid / Enteric fever": {
        "setting": "Community / Hospital",
        "likely_pathogens": ["Salmonella typhi","Salmonella paratyphi"],
        "first_line": [("Ciprofloxacin","500mg PO BD × 10-14d [UCG first-line]","W"),("Ceftriaxone","2g IV OD × 10-14d (severe/complicated)","W")],
        "second_line": [("Azithromycin","20mg/kg/d PO × 7d (WHO/international guideline preferred if FQ resistance suspected)","W")],
        "avoid": "Empiric ciprofloxacin (FQ DCS/resistance >40% in E. Africa). Empiric ampicillin/chloramphenicol/TMP-SMX (MDR rates high).",
        "notes": "Blood culture is gold standard (sensitivity ~60%). UCG 2023 recommends ciprofloxacin first-line. However, WHO and international guidelines increasingly prefer azithromycin due to rising FQ resistance and XDR S. Typhi (H58 lineage). In Uganda, ciprofloxacin remains acceptable per UCG pending local resistance data. C&S recommended. TCV vaccination recommended."
    },
    "Tuberculous meningitis": {
        "setting": "Hospital",
        "likely_pathogens": ["Mycobacterium tuberculosis"],
        "first_line": [("RHZE","(Rifampicin, INH, Pyrazinamide, Ethambutol) + Dexamethasone","A")],
        "second_line": [("If MDR-TB: Bedaquiline + Linezolid + Levofloxacin + Cycloserine","Specialist MDR regimen","R")],
        "avoid": "Ethambutol in children where visual monitoring impossible (use 3-drug RHZ). Delay in steroids.",
        "notes": "Dexamethasone reduces mortality (Thwaites et al. NEJM 2004). 12 months treatment (2 RHZE + 10 RH). GeneXpert on CSF (sensitivity ~60-80%). High mortality (~25-50%)."
    },
    "Cryptococcal meningitis (HIV)": {
        "setting": "Hospital",
        "likely_pathogens": ["Cryptococcus neoformans"],
        "first_line": [("Amphotericin B deoxycholate + Flucytosine","1 week induction, then Fluconazole 800mg/d consolidation","R/W")],
        "second_line": [("Amphotericin B + Fluconazole","(if flucytosine unavailable)","R/W"),("Fluconazole 1200mg/d + Flucytosine","(if amphotericin unavailable)","W")],
        "avoid": "Fluconazole monotherapy for induction (inferior outcomes). Immediate ART initiation (wait 4-6 weeks — COAT trial).",
        "notes": "WHO 2022 guidelines. CrAg screening at CD4<200 prevents disease. Therapeutic LP (opening pressure management) critical. ART delay 4-6 weeks post-diagnosis."
    },
    # ── ADDED SYNDROMES (HIGH PRIORITY) ──────────────────────────────────
    "Acute bacterial diarrhoea / dysentery": {
        "setting": "Community / Hospital",
        "likely_pathogens": ["Shigella spp.","Non-typhoidal Salmonella (iNTS)","Campylobacter jejuni","Vibrio cholerae","Enterotoxigenic E. coli (ETEC)"],
        "first_line": [("ORS + Zinc","(all cases — rehydration is primary treatment)","A"),("Ciprofloxacin","500mg PO BD × 3d (bloody diarrhoea/dysentery)","W")],
        "second_line": [("Azithromycin","500mg PO OD × 3d (if FQ resistance suspected/children)","W"),("Ceftriaxone","1-2g IV OD (severe/systemic signs)","W")],
        "avoid": "Antibiotics for watery diarrhoea without bloody stool or systemic toxicity (most cases viral/self-limiting). Loperamide in dysentery (risk of toxic megacolon).",
        "notes": "ORS is the primary intervention for ALL acute diarrhoea. Antibiotics only for: bloody diarrhoea (dysentery), cholera (suspected), traveller's diarrhoea (moderate-severe), or immunocompromised. Stool culture if possible. Shigella FQ resistance rising in E. Africa."
    },
    "Pelvic inflammatory disease (PID)": {
        "setting": "Community / Hospital",
        "likely_pathogens": ["Neisseria gonorrhoeae","Chlamydia trachomatis","Anaerobes (Bacteroides, Prevotella)","Mycoplasma genitalium","Enterobacteriaceae"],
        "first_line": [("Ceftriaxone 500mg IM stat + Doxycycline 100mg PO BD × 14d + Metronidazole 400mg PO BD × 14d","(WHO/CDC standard triple therapy)","W+A+A")],
        "second_line": [("Ceftriaxone + Azithromycin 1g PO stat + Metronidazole","(if doxycycline unavailable/not tolerated)","W+W+A")],
        "avoid": "Fluoroquinolone monotherapy (poor Chlamydia activity, rising gonococcal resistance). Delay in treatment (risk of tubal damage and infertility).",
        "notes": "Low threshold for treatment — clinical diagnosis (lower abdominal pain + cervical motion tenderness) sufficient to start empiric therapy. IUD removal NOT required if clinically improving. Partner treatment essential. Test for HIV."
    },
    "Diabetic foot infection": {
        "setting": "Community / Hospital",
        "likely_pathogens": ["Staphylococcus aureus","Group A Streptococcus","Enterobacteriaceae (E. coli, Proteus, Klebsiella)","Anaerobes (Bacteroides)","Pseudomonas aeruginosa (chronic wounds)"],
        "first_line": [("Amoxicillin-clavulanate","625mg PO TDS (mild, no MRSA risk)","A"),("Ceftriaxone + Metronidazole","IV (moderate-severe)","W+A")],
        "second_line": [("Piperacillin-tazobactam","4.5g IV Q6h (severe/limb-threatening)","W"),("Meropenem + Vancomycin","(severe with MRSA + ESBL risk)","W+W")],
        "avoid": "Antibiotics alone without wound assessment, debridement, and offloading. Topical antibiotics on deep infections.",
        "notes": "IDSA/IWGDF classification: mild (superficial, <2cm cellulitis), moderate (deeper/larger), severe (systemic toxicity/sepsis). Probe-to-bone test for osteomyelitis. X-ray foot. Vascular assessment (ABI). Multidisciplinary approach: surgery, medicine, podiatry. Duration: mild 1-2 weeks, moderate 2-3 weeks, with osteomyelitis 4-6 weeks."
    },
    "Septic arthritis": {
        "setting": "Hospital",
        "likely_pathogens": ["Staphylococcus aureus","Streptococcus spp.","Neisseria gonorrhoeae (young adults)","Gram-negative rods (elderly, immunocompromised)","Salmonella spp. (sickle cell disease)"],
        "first_line": [("Flucloxacillin","2g IV Q6h","A"),("+ Ceftriaxone","if Gram-negative or gonococcal suspected","W")],
        "second_line": [("Vancomycin","IV (if MRSA risk)","W"),("Clindamycin","IV (penicillin allergy)","A")],
        "avoid": "Delay in joint aspiration and washout. Empiric treatment without synovial fluid analysis (Gram stain, culture, cell count, crystals).",
        "notes": "EMERGENCY — joint destruction within 24-48h if untreated. Aspirate BEFORE antibiotics. Synovial WBC >50,000 with >90% neutrophils = septic until proven otherwise. Surgical washout/arthroscopy for hip and most joints. In sickle cell disease, consider Salmonella (ceftriaxone covers). Gonococcal arthritis: ceftriaxone 1g IV OD, often migratory polyarthralgia."
    },
    "Osteomyelitis (acute)": {
        "setting": "Hospital",
        "likely_pathogens": ["Staphylococcus aureus","Streptococcus spp.","Enterobacteriaceae","Salmonella spp. (sickle cell)","Pseudomonas aeruginosa (puncture wounds)"],
        "first_line": [("Flucloxacillin","2g IV Q6h × 2 weeks then PO × 4 weeks","A"),("Ceftriaxone","2g IV OD (if Gram-negative/Salmonella suspected)","W")],
        "second_line": [("Vancomycin","IV (MRSA)","W"),("Clindamycin","600mg IV/PO (penicillin allergy, good bone penetration)","A")],
        "avoid": "Premature switch to oral without clinical improvement and falling CRP. Inadequate duration (<4 weeks total).",
        "notes": "Blood cultures (positive in ~50%). MRI is gold standard imaging (>90% sensitivity). Bone biopsy for culture if blood cultures negative. Sickle cell: Salmonella most common (not Staph) — use ceftriaxone. Duration: 4-6 weeks total (can switch IV→PO after 1-2 weeks if improving, CRP falling — OVIVA trial)."
    },
    "Acute otitis media": {
        "setting": "Community",
        "likely_pathogens": ["Streptococcus pneumoniae","Haemophilus influenzae","Moraxella catarrhalis","Group A Streptococcus"],
        "first_line": [("Amoxicillin","80-90mg/kg/d PO (children) or 500mg TDS (adults) × 5-7d","A")],
        "second_line": [("Amoxicillin-clavulanate","(if no response at 48-72h — beta-lactamase producing H. influenzae)","A"),("Ceftriaxone","50mg/kg IM × 3d (vomiting child)","W")],
        "avoid": "Immediate antibiotics for mild AOM in children >2 years with unilateral disease (watchful waiting 48-72h acceptable per AAP). Ear drops alone for AOM (unless perforated TM with otorrhoea).",
        "notes": "Most common reason for antibiotic prescribing in children. High-dose amoxicillin overcomes pneumococcal intermediate resistance. Tympanocentesis if: treatment failure, severe/recurrent, immunocompromised, neonate. Prevent with PCV13 vaccination."
    },
    "Acute pharyngitis / tonsillitis": {
        "setting": "Community",
        "likely_pathogens": ["Group A Streptococcus (S. pyogenes)","Viral (majority — EBV, adenovirus, rhinovirus)","Fusobacterium necrophorum (adolescents)"],
        "first_line": [("Phenoxymethylpenicillin (Pen V)","500mg PO BD × 10d","A"),("Amoxicillin","50mg/kg OD × 10d (children — better taste)","A")],
        "second_line": [("Azithromycin","500mg D1 then 250mg D2-5 (penicillin allergy)","W"),("Cefalexin","500mg PO BD × 10d (non-anaphylactic penicillin allergy)","A")],
        "avoid": "Amoxicillin if EBV/glandular fever suspected (causes widespread maculopapular rash). Antibiotics for confirmed viral pharyngitis.",
        "notes": "CRITICAL in East Africa: untreated GAS pharyngitis → rheumatic fever → rheumatic heart disease. Centor/McIsaac score guides need for antibiotics (≥3 = treat or test). Full 10-day penicillin course required for rheumatic fever prevention (shorter courses do NOT eradicate GAS). If recurrent, consider tonsillectomy referral."
    },
    "Catheter-associated UTI (CAUTI)": {
        "setting": "Hospital",
        "likely_pathogens": ["Escherichia coli","Klebsiella pneumoniae","Pseudomonas aeruginosa","Enterococcus spp.","Candida spp. (prolonged catheterisation)"],
        "first_line": [("Remove/replace catheter","(ESSENTIAL first step)","—"),("Ceftriaxone","1-2g IV OD (pending culture)","W")],
        "second_line": [("Piperacillin-tazobactam","4.5g IV Q6h (Pseudomonas risk)","W"),("Meropenem","1g IV Q8h (ESBL risk/prior resistant cultures)","W")],
        "avoid": "Treating asymptomatic bacteriuria in catheterised patients (except pre-urological procedures or pregnancy). Catheter change alone without systemic antibiotics if symptomatic.",
        "notes": "Remove or replace catheter BEFORE starting antibiotics — biofilm on old catheter harbours resistant organisms. Duration: 7 days if prompt response, 10-14 days if delayed. Candiduria: remove catheter; fluconazole only if symptomatic/high-risk. Prevention: remove catheter as early as possible (nurse-led reminder protocols reduce CAUTI 30-50%)."
    },
    "Surgical site infection (SSI)": {
        "setting": "Hospital",
        "likely_pathogens": ["Staphylococcus aureus","Coagulase-negative Staphylococcus","Escherichia coli","Enterococcus spp.","Anaerobes (abdominal surgery)","Pseudomonas aeruginosa (prolonged hospitalisation)"],
        "first_line": [("Flucloxacillin","500mg-1g PO/IV (superficial, MSSA)","A"),("Amoxicillin-clavulanate","1.2g IV TDS (deep incisional, abdominal surgery)","A")],
        "second_line": [("Ceftriaxone + Metronidazole","IV (deep/organ-space, abdominal)","W+A"),("Vancomycin","IV (MRSA risk or prosthetic material)","W")],
        "avoid": "Prolonged prophylactic antibiotics postoperatively (single dose ± 24h is sufficient for prophylaxis; extending beyond does NOT prevent SSI and drives resistance).",
        "notes": "Classification: superficial incisional (skin/subcutaneous), deep incisional (fascia/muscle), organ/space. Wound culture before antibiotics. Open and drain if collection present. SSI prevention: appropriate prophylaxis timing (within 60 min of incision), skin prep, normothermia, glycaemic control."
    },
    "Puerperal sepsis / endometritis": {
        "setting": "Hospital (Obstetric)",
        "likely_pathogens": ["Group A Streptococcus","Escherichia coli","Anaerobes (Bacteroides, Peptostreptococcus)","Group B Streptococcus","Enterococcus spp.","Staphylococcus aureus (wound infection)"],
        "first_line": [("Ampicillin + Gentamicin + Metronidazole","IV triple therapy (WHO recommended)","A+A+A")],
        "second_line": [("Amoxicillin-clavulanate","1.2g IV TDS (if triple therapy unavailable)","A"),("Piperacillin-tazobactam","4.5g IV Q6h (severe/no response)","W")],
        "avoid": "Delay — puerperal GAS sepsis can kill within hours. Oral antibiotics alone for febrile postpartum patients with signs of sepsis.",
        "notes": "OBSTETRIC EMERGENCY. WHO defines puerperal sepsis as infection of the genital tract occurring between rupture of membranes/labour and 42 days postpartum. High vaginal swab + blood cultures before antibiotics. Examine for retained products (ultrasound). If no improvement at 48h: CT abdomen/pelvis for pelvic abscess. Thromboprophylaxis — sepsis is a VTE risk factor."
    },
    # ── ADDED SYNDROMES (MEDIUM PRIORITY) ────────────────────────────────
    "Acute cholecystitis / cholangitis": {
        "setting": "Hospital (Surgical)",
        "likely_pathogens": ["Escherichia coli","Klebsiella pneumoniae","Enterococcus spp.","Anaerobes (Bacteroides, Clostridium)"],
        "first_line": [("Ceftriaxone + Metronidazole","IV (cholecystitis)","W+A"),("Piperacillin-tazobactam","4.5g IV Q6h (cholangitis — broader cover)","W")],
        "second_line": [("Meropenem","1g IV Q8h (severe/ESBL risk)","W"),("Amoxicillin-clavulanate","1.2g IV TDS (mild cholecystitis)","A")],
        "avoid": "Conservative management alone in acute cholangitis with obstruction (requires urgent ERCP/drainage).",
        "notes": "Tokyo Guidelines severity grading. Cholecystitis: early cholecystectomy (within 72h) preferred. Cholangitis (Charcot's triad: fever, jaundice, RUQ pain; Reynolds' pentad adds confusion + hypotension): urgent biliary drainage (ERCP) within 24h. Blood cultures before antibiotics."
    },
    "Pyogenic liver abscess": {
        "setting": "Hospital",
        "likely_pathogens": ["Escherichia coli","Klebsiella pneumoniae (hypervirulent K1/K2)","Anaerobes","Streptococcus milleri group","Entamoeba histolytica (amoebic — distinguish)"],
        "first_line": [("Ceftriaxone + Metronidazole","IV (covers bacterial + amoebic empirically)","W+A")],
        "second_line": [("Piperacillin-tazobactam","4.5g IV Q6h","W"),("Meropenem","(if ESBL risk)","W")],
        "avoid": "Assuming amoebic without investigation — pyogenic abscesses require drainage; amoebic usually respond to metronidazole alone.",
        "notes": "Distinguish amoebic vs pyogenic: amoebic = single, right lobe, Entamoeba serology+; pyogenic = may be multiple, often with biliary source. Percutaneous drainage for abscesses >5cm or not responding to antibiotics at 72h. Duration: 4-6 weeks total (IV then PO step-down). Blood cultures positive in 50%."
    },
    "Brain abscess": {
        "setting": "Hospital (Neurosurgical)",
        "likely_pathogens": ["Streptococcus milleri group","Staphylococcus aureus (trauma/post-surgical)","Anaerobes","Enterobacteriaceae","Nocardia (immunocompromised)","Toxoplasma gondii (HIV — ring-enhancing, treat empirically)"],
        "first_line": [("Ceftriaxone 2g IV BD + Metronidazole 500mg IV TDS","(empiric, covers strep + anaerobes)","W+A")],
        "second_line": [("+ Vancomycin","(if post-surgical or MRSA risk)","W"),("+ TMP-SMX","(if Nocardia suspected — immunocompromised)","A")],
        "avoid": "Lumbar puncture (risk of herniation if mass effect). Steroids before diagnosis confirmed (may shrink lymphoma/toxoplasmosis and confound diagnosis).",
        "notes": "CT/MRI with contrast: ring-enhancing lesion(s). Source: contiguous (sinusitis, otitis, dental) or haematogenous (endocarditis, lung abscess). Stereotactic aspiration for diagnosis + drainage if >2.5cm. Duration: 6-8 weeks IV antibiotics. In HIV: treat empirically for toxoplasmosis first (pyrimethamine + sulfadiazine) — biopsy if no response at 2 weeks."
    },
    "Late-onset neonatal sepsis": {
        "setting": "Hospital (Neonatal unit)",
        "likely_pathogens": ["Coagulase-negative Staphylococcus (CLABSI)","Staphylococcus aureus","Klebsiella pneumoniae (ESBL — major E. Africa)","Escherichia coli","Pseudomonas aeruginosa","Candida spp. (VLBW, prolonged antibiotics)"],
        "first_line": [("Vancomycin + Gentamicin","IV (covers CoNS/MRSA + Gram-negatives)","W+A")],
        "second_line": [("Meropenem","IV (if ESBL Klebsiella suspected — common E. African NICUs)","W"),("+ Fluconazole","(if Candida risk: VLBW, prolonged broad-spectrum antibiotics, CVC)","W")],
        "avoid": "Empiric ceftriaxone in neonates receiving calcium infusions (fatal ceftriaxone-calcium precipitate). Prolonged empiric antibiotics without positive cultures (drives NEC and fungal infection).",
        "notes": "Onset >72h of life. Different pathogen spectrum from early-onset: nosocomial organisms, line-related. ESBL Klebsiella is the leading killer in many E. African NICUs. CVC removal if line sepsis suspected. Blood culture before antibiotics. Duration: 7-10 days for bacteraemia; 14-21 days for meningitis."
    },
    "Bite wound infection (human / animal)": {
        "setting": "Community / Hospital",
        "likely_pathogens": ["Pasteurella multocida (cat/dog bite)","Eikenella corrodens (human bite)","Staphylococcus aureus","Streptococcus spp.","Anaerobes (Fusobacterium, Bacteroides)","Capnocytophaga canimorsus (dog bite, asplenic patients)"],
        "first_line": [("Amoxicillin-clavulanate","625mg PO TDS × 5-7d (prophylaxis or established infection)","A")],
        "second_line": [("Doxycycline + Metronidazole","(penicillin allergy)","A+A"),("Ceftriaxone + Metronidazole","IV (severe/systemic)","W+A")],
        "avoid": "Primary closure of bite wounds (except face — cosmetic closure with prophylactic antibiotics). Flucloxacillin monotherapy (misses Pasteurella and anaerobes).",
        "notes": "Irrigate copiously with saline. Tetanus prophylaxis. Rabies PEP for all animal bites in endemic areas (Uganda is endemic — PEP indicated for ALL dog/cat bites unless animal observed healthy 10 days). Human bites: high infection risk, X-ray for tooth fragments, assess for clenched-fist injury (fight bite). Prophylactic antibiotics for: cat bites (high infection rate), hand bites, immunocompromised, deep puncture wounds."
    },
    "Necrotising fasciitis": {
        "setting": "Hospital (Surgical EMERGENCY)",
        "likely_pathogens": ["Group A Streptococcus (Type II — monomicrobial)","Mixed aerobic-anaerobic (Type I — polymicrobial, often post-surgical/diabetic)","Clostridium perfringens (gas gangrene)","Vibrio vulnificus (saltwater exposure)"],
        "first_line": [("Meropenem 1g IV Q8h + Clindamycin 600mg IV Q8h + Vancomycin","(empiric triple — covers polymicrobial + GAS + MRSA + toxin suppression)","W+A+W")],
        "second_line": [("Piperacillin-tazobactam + Clindamycin + Vancomycin","(alternative broad-spectrum)","W+A+W")],
        "avoid": "ANY delay in surgical exploration. Waiting for imaging confirmation (clinical diagnosis sufficient — pain out of proportion, crepitus, rapid spread, systemic toxicity). Antibiotics alone without surgery.",
        "notes": "SURGICAL EMERGENCY — mortality 20-40% even with treatment. LRINEC score aids diagnosis (WBC, Hb, Na, glucose, creatinine, CRP). Key signs: pain out of proportion, woody induration, skin necrosis, crepitus, haemodynamic instability. Immediate surgical exploration and radical debridement. Clindamycin suppresses toxin production (Eagle effect). Often requires repeat debridement every 24-48h. IVIG may benefit in streptococcal toxic shock."
    },
    "Lung abscess / empyema": {
        "setting": "Hospital",
        "likely_pathogens": ["Anaerobes (Fusobacterium, Peptostreptococcus, Bacteroides)","Streptococcus milleri group","Staphylococcus aureus","Klebsiella pneumoniae (alcoholism, diabetes)"],
        "first_line": [("Amoxicillin-clavulanate","1.2g IV TDS","A"),("Clindamycin","600mg IV Q8h (alternative — good anaerobic + lung penetration)","A")],
        "second_line": [("Ceftriaxone + Metronidazole","IV","W+A"),("Meropenem","(if ESBL/hospital-acquired)","W")],
        "avoid": "Percutaneous drainage of lung abscess (risk of bronchopleural fistula and empyema — most resolve with antibiotics alone). CT-guided drainage only if peripheral and not communicating with bronchus.",
        "notes": "Lung abscess: most respond to 4-6 weeks antibiotics alone (postural drainage, chest physiotherapy). Empyema: requires chest tube drainage + antibiotics. BTS guidelines: pH <7.2, positive Gram stain/culture, or frank pus = tube drainage. If loculated, consider intrapleural fibrinolytics (alteplase + DNase) or VATS. Exclude underlying malignancy (bronchoscopy if no risk factors for aspiration)."
    },
    "Catheter-related bloodstream infection (CRBSI)": {
        "setting": "Hospital / ICU",
        "likely_pathogens": ["Coagulase-negative Staphylococcus","Staphylococcus aureus","Enterococcus spp.","Candida spp.","Gram-negative rods (Klebsiella, Pseudomonas, Acinetobacter)"],
        "first_line": [("Remove CVC","(ESSENTIAL unless tunnelled/no alternative access)","—"),("Vancomycin","IV (covers CoNS + MRSA)","W")],
        "second_line": [("+ Cefepime or Meropenem","(add if Gram-negative suspected or ICU patient)","W"),("+ Caspofungin/Fluconazole","(if Candida risk — TPN, prolonged antibiotics, prior colonisation)","W")],
        "avoid": "Retaining catheter for S. aureus or Candida CRBSI (catheter MUST be removed). Antibiotic lock therapy alone without systemic antibiotics.",
        "notes": "Differential time to positivity (DTP): blood from CVC positive ≥2h before peripheral = CRBSI. IDSA guidelines: remove catheter for S. aureus, Candida, P. aeruginosa, mycobacteria. S. aureus CRBSI: minimum 4 weeks IV antibiotics + echocardiogram (endocarditis in ~25%). CoNS: may salvage catheter with antibiotic lock if uncomplicated."
    },
    "Spontaneous bacterial peritonitis (SBP)": {
        "setting": "Hospital",
        "likely_pathogens": ["Escherichia coli","Klebsiella pneumoniae","Streptococcus pneumoniae","Enterococcus spp."],
        "first_line": [("Ceftriaxone","2g IV OD","W"),("+ IV Albumin","1.5g/kg day 1, 1g/kg day 3 (hepatorenal syndrome prevention)","—")],
        "second_line": [("Ciprofloxacin","400mg IV BD (if cephalosporin allergy)","W"),("Meropenem","(if healthcare-associated or FQ prophylaxis failure)","W")],
        "avoid": "Aminoglycosides (nephrotoxic in cirrhosis — contraindicated). Delay in diagnostic paracentesis.",
        "notes": "Diagnostic paracentesis for ALL cirrhotic patients with ascites admitted to hospital. SBP = ascitic fluid PMN ≥250/mm³. Start antibiotics immediately after paracentesis. Duration: 5 days. Secondary prophylaxis: norfloxacin 400mg OD or co-trimoxazole. Albumin reduces renal impairment and mortality (Sort et al. NEJM 1999). Distinguish from secondary peritonitis (multiple organisms, very high PMN, protein >1g/dL — requires surgical evaluation)."
    },
    "Febrile neutropenia": {
        "setting": "Hospital (Oncology / Haematology)",
        "likely_pathogens": ["Escherichia coli","Klebsiella pneumoniae","Pseudomonas aeruginosa","Staphylococcus aureus","Coagulase-negative Staphylococcus","Streptococcus viridans"],
        "first_line": [("Piperacillin-tazobactam","4.5g IV Q6h (IDSA first-line)","W"),("Meropenem","1g IV Q8h (if ESBL colonisation/high local ESBL rates)","W")],
        "second_line": [("+ Vancomycin","(add if: MRSA risk, line infection, mucositis, skin/soft tissue infection, haemodynamic instability)","W"),("+ Caspofungin","(add antifungal if fever persists >4-7 days on broad-spectrum antibiotics)","W")],
        "avoid": "Delay — antibiotics within 1 hour of fever onset (mortality doubles with each hour delay). Oral empiric therapy only in MASCC score ≥21 (low risk) with close follow-up.",
        "notes": "Defined as: ANC <500/mm³ (or expected to fall) + single temperature ≥38.3°C or sustained ≥38.0°C over 1 hour. MASCC score for risk stratification. Blood cultures (2 sets: one from CVC if present, one peripheral) BEFORE antibiotics. Duration: until afebrile ≥48h AND ANC recovering (>500). G-CSF not routinely recommended for treatment but consider if high-risk features."
    },
    "Acute prostatitis": {
        "setting": "Community / Hospital",
        "likely_pathogens": ["Escherichia coli","Klebsiella pneumoniae","Proteus mirabilis","Pseudomonas aeruginosa","Enterococcus spp.","Neisseria gonorrhoeae / Chlamydia (young men — STI-associated)"],
        "first_line": [("Ciprofloxacin","500mg PO BD × 4 weeks (good prostatic penetration)","W"),("Ceftriaxone","2g IV OD (if systemically unwell, then switch to cipro when improving)","W")],
        "second_line": [("Co-trimoxazole","960mg PO BD × 4 weeks (if susceptible)","A"),("Doxycycline","100mg PO BD × 4 weeks (if STI-associated)","A")],
        "avoid": "Short courses (<2 weeks) — risk of chronic prostatitis/relapse. Nitrofurantoin (no prostatic penetration).",
        "notes": "Acutely swollen, exquisitely tender prostate on DRE (gentle exam only — do NOT massage). Urine culture + STI screen in young men. If urinary retention: suprapubic catheter (urethral catheter causes pain/risk). 4-week course minimum required for prostatic penetration. If abscess suspected (failure to improve): transrectal ultrasound, drainage if >2cm."
    },
    "Syphilis": {
        "setting": "Community / STI clinic",
        "likely_pathogens": ["Treponema pallidum"],
        "first_line": [("Benzathine penicillin G","2.4 MU IM single dose (primary/secondary)","A"),("Benzathine penicillin G","2.4 MU IM weekly × 3 weeks (tertiary/latent >1yr)","A")],
        "second_line": [("Doxycycline","100mg PO BD × 14d (penicillin allergy, non-pregnant)","A"),("Ceftriaxone","1g IM/IV OD × 10-14d (alternative)","W")],
        "avoid": "Azithromycin (rising resistance 30-90% in some populations — NOT recommended for syphilis). Oral penicillin (inadequate levels).",
        "notes": "No documented penicillin resistance in T. pallidum — ever. Benzathine penicillin G is the ONLY option in pregnancy (desensitise if allergic). Neurosyphilis: benzylpenicillin 3-4 MU IV Q4h × 14d. Screen all pregnant women (RPR/VDRL). Jarisch-Herxheimer reaction may occur in first 24h of treatment. UCG 2023 aligned: benzathine penicillin G first-line."
    },
    "Peptic ulcer disease (H. pylori)": {
        "setting": "Community / Hospital",
        "likely_pathogens": ["Helicobacter pylori"],
        "first_line": [("PPI + Amoxicillin 1g BD + Clarithromycin 500mg BD","14 days (standard triple therapy)","A+W")],
        "second_line": [("PPI + Bismuth subsalicylate + Metronidazole 500mg TDS + Tetracycline 500mg QDS","14 days (bismuth quadruple — if clarithromycin resistance suspected)","A"),("PPI + Amoxicillin + Metronidazole + Clarithromycin","14 days (concomitant quadruple)","A+W")],
        "avoid": "Clarithromycin-based triple therapy if local clarithromycin resistance >15%. Metronidazole monotherapy (high resistance).",
        "notes": "Test for H. pylori before treating (stool antigen or urea breath test). Confirm eradication 4 weeks after treatment. In Africa, metronidazole resistance is very high (>60%) — consider bismuth quadruple or concomitant therapy. UCG 2023: triple therapy with PPI + amoxicillin + clarithromycin."
    },
    "Cholera": {
        "setting": "Community / Cholera treatment centre",
        "likely_pathogens": ["Vibrio cholerae"],
        "first_line": [("ORS/IV Ringer lactate","Rehydration is the PRIMARY treatment","A"),("Doxycycline","300mg single dose (adults) OR 100mg BD × 3d","A")],
        "second_line": [("Azithromycin","1g single dose (children, pregnant women)","W"),("Ciprofloxacin","1g single dose (alternative)","W")],
        "avoid": "Antibiotics WITHOUT adequate rehydration (rehydration saves lives, antibiotics are adjunctive). TMP-SMX (high resistance).",
        "notes": "NOTIFIABLE DISEASE. ORS/IV fluids are life-saving — antibiotics reduce duration by ~50% and shedding. Zinc supplementation in children. Oral cholera vaccine for outbreak control. UCG 2023 aligned."
    },
    "Pertussis (whooping cough)": {
        "setting": "Community",
        "likely_pathogens": ["Bordetella pertussis"],
        "first_line": [("Azithromycin","10mg/kg/d × 5d (infants); 500mg D1 then 250mg D2-5 (adults/children)","W")],
        "second_line": [("Erythromycin","500mg QDS × 14d (adults)","A"),("Co-trimoxazole","(if macrolide contraindicated)","A")],
        "avoid": "Starting antibiotics >3 weeks after cough onset (no clinical benefit, but still reduces transmission). Cephalosporins (not effective against Bordetella).",
        "notes": "Antibiotics reduce transmission but do NOT shorten clinical course if started late. Most effective in catarrhal phase (first 1-2 weeks). Infants <6 months: hospitalise (apnoea risk). Prophylaxis for close contacts (especially pregnant women near term and infants). DPT vaccination is primary prevention."
    },
}
# ── 1D. ICD-10 BACTERIAL INFECTION CODE DATABASE ─────────────────────────
# Maps ICD-10-CM codes → clinical context, likely pathogens, and linked syndrome
# Source: WHO ICD-10 2019 edition; clinical pathogen mapping from MANDELL 9th Ed
ICD_10 = {
    # ── INTESTINAL INFECTIONS ──
    "A00":"Cholera|Intestinal|Vibrio cholerae|ORS + doxycycline/azithromycin; fluid resuscitation",
    "A01.0":"Typhoid fever|Intestinal|Salmonella typhi|Typhoid / Enteric fever",
    "A01.1":"Paratyphoid fever A|Intestinal|Salmonella paratyphi A|Typhoid / Enteric fever",
    "A01.4":"Enteric fever, unspecified|Intestinal|Salmonella typhi, S. paratyphi|Typhoid / Enteric fever",
    "A02.0":"Salmonella enteritis|Intestinal|Non-typhoidal Salmonella|Usually self-limiting; antibiotics if invasive/immunocompromised",
    "A02.1":"Salmonella sepsis|Bloodstream|Non-typhoidal Salmonella (iNTS)|Sepsis / Septic shock (undifferentiated)",
    "A03":"Shigellosis|Intestinal|Shigella spp.|Ciprofloxacin or azithromycin; rising FQ resistance",
    "A04.5":"Campylobacter enteritis|Intestinal|Campylobacter jejuni|Usually self-limiting; azithromycin if severe",
    "A04.7":"Clostridioides difficile colitis|Intestinal|Clostridioides difficile|Oral vancomycin first-line (IDSA 2021)",
    "A04.8":"Other bacterial intestinal infections|Intestinal|Various|Culture-directed",
    "A09":"Infectious gastroenteritis/colitis|Intestinal|Multiple bacterial/viral|ORS; antibiotics only if bacterial confirmed/suspected",
    # ── TUBERCULOSIS ──
    "A15.0":"TB of lung, confirmed by sputum microscopy|Respiratory|Mycobacterium tuberculosis|RHZE standard regimen",
    "A15.3":"TB of lung, confirmed by culture|Respiratory|Mycobacterium tuberculosis|RHZE; DST-guided if MDR suspected",
    "A15.7":"Primary respiratory TB, bacteriologically confirmed|Respiratory|Mycobacterium tuberculosis|RHZE standard",
    "A16.0":"TB of lung, bacteriology negative|Respiratory|Mycobacterium tuberculosis|Clinical TB; RHZE",
    "A16.2":"TB of lung, without mention of confirmation|Respiratory|Mycobacterium tuberculosis|RHZE standard",
    "A17.0":"Tuberculous meningitis|CNS|Mycobacterium tuberculosis|Tuberculous meningitis",
    "A17.1":"Meningeal tuberculoma|CNS|Mycobacterium tuberculosis|Tuberculous meningitis",
    "A18.0":"TB of bones and joints|Musculoskeletal|Mycobacterium tuberculosis|RHZE 12 months; surgical if spinal cord compromise",
    "A18.1":"TB of genitourinary system|Genitourinary|Mycobacterium tuberculosis|RHZE standard",
    "A18.2":"TB peripheral lymphadenopathy|Lymphatic|Mycobacterium tuberculosis|RHZE standard 6 months",
    "A18.3":"TB of intestines/peritoneum|Abdominal|Mycobacterium tuberculosis|RHZE; exclude surgical emergency",
    "A18.4":"TB of skin and subcutaneous tissue|Skin|Mycobacterium tuberculosis|RHZE standard",
    "A18.5":"TB of eye|Ophthalmic|Mycobacterium tuberculosis|RHZE + ophthalmology input",
    "A19.0":"Acute miliary TB, single specified site|Disseminated|Mycobacterium tuberculosis|RHZE; consider TB meningitis workup",
    "A19.1":"Acute miliary TB, multiple sites|Disseminated|Mycobacterium tuberculosis|RHZE; high mortality",
    # ── ZOONOTIC / SPECIFIC BACTERIA ──
    "A22.0":"Cutaneous anthrax|Skin|Bacillus anthracis|Ciprofloxacin or doxycycline",
    "A23.0":"Brucellosis due to B. melitensis|Systemic|Brucella melitensis|Doxycycline + streptomycin/gentamicin 6 weeks",
    "A23.9":"Brucellosis, unspecified|Systemic|Brucella spp.|Doxycycline + rifampicin or aminoglycoside",
    "A27.0":"Leptospirosis (Weil disease)|Systemic|Leptospira interrogans|Penicillin G or ceftriaxone IV (severe); doxycycline (mild)",
    "A28.2":"Extraintestinal yersiniosis|Systemic|Yersinia enterocolitica|Ciprofloxacin or TMP-SMX + aminoglycoside",
    "A32.0":"Cutaneous listeriosis|Skin|Listeria monocytogenes|Ampicillin; add gentamicin if severe",
    "A32.1":"Listerial meningitis|CNS|Listeria monocytogenes|Bacterial meningitis (community-acquired)",
    "A32.7":"Listerial sepsis|Bloodstream|Listeria monocytogenes|Ampicillin + gentamicin IV",
    # ── STREPTOCOCCAL / STAPHYLOCOCCAL ──
    "A38":"Scarlet fever|Systemic|Group A Streptococcus|Penicillin V or amoxicillin PO",
    "A39.0":"Meningococcal meningitis|CNS|Neisseria meningitidis|Bacterial meningitis (community-acquired)",
    "A39.1":"Waterhouse-Friderichsen syndrome|Bloodstream|Neisseria meningitidis|Ceftriaxone IV + intensive care",
    "A39.2":"Acute meningococcaemia|Bloodstream|Neisseria meningitidis|Ceftriaxone IV emergency",
    "A39.4":"Meningococcaemia, unspecified|Bloodstream|Neisseria meningitidis|Ceftriaxone IV",
    "A40.0":"Sepsis due to Group A Streptococcus|Bloodstream|Group A Streptococcus (S. pyogenes)|Penicillin + clindamycin (toxin suppression)",
    "A40.1":"Sepsis due to Group B Streptococcus|Bloodstream|Group B Streptococcus (S. agalactiae)|Penicillin G or ampicillin IV",
    "A40.3":"Sepsis due to S. pneumoniae|Bloodstream|Streptococcus pneumoniae|Ceftriaxone IV",
    "A41.0":"Sepsis due to S. aureus|Bloodstream|Staphylococcus aureus|Sepsis / Septic shock (undifferentiated)",
    "A41.1":"Sepsis due to other Staphylococcus|Bloodstream|Coagulase-negative Staphylococcus|Vancomycin IV; remove line if CLABSI",
    "A41.5":"Sepsis due to other Gram-negative organisms|Bloodstream|Escherichia coli, Klebsiella pneumoniae|Sepsis / Septic shock (undifferentiated)",
    "A41.9":"Sepsis, unspecified|Bloodstream|Unknown|Sepsis / Septic shock (undifferentiated)",
    "A46":"Erysipelas|Skin|Group A Streptococcus|Skin and soft tissue infection (SSTI)",
    "A48.0":"Gas gangrene|Skin|Clostridium perfringens|Penicillin + clindamycin + surgical debridement; EMERGENCY",
    "A48.1":"Legionnaires disease|Respiratory|Legionella pneumophila|Azithromycin or levofloxacin",
    "A49.0":"Staphylococcal infection, unspecified|Various|Staphylococcus aureus|Skin and soft tissue infection (SSTI)",
    "A49.1":"Streptococcal infection, unspecified|Various|Streptococcus spp.|Culture-directed; penicillin first-line",
    "A49.9":"Bacterial infection, unspecified|Various|Unknown|Culture before antibiotics; empiric per site of infection",
    # ── STIs ──
    "A50.0":"Early congenital syphilis, symptomatic|Neonatal|Treponema pallidum|Syphilis",
    "A50.9":"Congenital syphilis, unspecified|Neonatal|Treponema pallidum|Syphilis",
    "A51.0":"Primary genital syphilis (chancre)|Genitourinary|Treponema pallidum|Syphilis",
    "A51.3":"Secondary syphilis of skin and mucous membranes|Systemic|Treponema pallidum|Syphilis",
    "A52.1":"Symptomatic neurosyphilis|CNS|Treponema pallidum|Syphilis",
    "A52.3":"Cardiovascular syphilis|Cardiovascular|Treponema pallidum|Syphilis",
    "A53.9":"Syphilis, unspecified|Systemic|Treponema pallidum|Syphilis",
    "A37.0":"Whooping cough due to B. pertussis|Respiratory|Bordetella pertussis|Pertussis (whooping cough)",
    "A37.9":"Whooping cough, unspecified|Respiratory|Bordetella pertussis|Pertussis (whooping cough)",
    "A00.0":"Cholera due to V. cholerae O1|Intestinal|Vibrio cholerae|Cholera",
    "A00.9":"Cholera, unspecified|Intestinal|Vibrio cholerae|Cholera",
    "K25.9":"Gastric ulcer (H. pylori associated)|Abdominal|Helicobacter pylori|Peptic ulcer disease (H. pylori)",
    "K26.9":"Duodenal ulcer (H. pylori associated)|Abdominal|Helicobacter pylori|Peptic ulcer disease (H. pylori)",
    "A54.0":"Gonococcal urethritis/cervicitis|Genitourinary|Neisseria gonorrhoeae|Sexually transmitted infection (urethral/cervical discharge)",
    "A54.1":"Gonococcal pelvic peritonitis/PID|Genitourinary|Neisseria gonorrhoeae|Ceftriaxone + doxycycline + metronidazole",
    "A54.2":"Gonococcal pelviperitonitis|Genitourinary|Neisseria gonorrhoeae|As PID regimen",
    "A54.9":"Gonococcal infection, unspecified|Genitourinary|Neisseria gonorrhoeae|Sexually transmitted infection (urethral/cervical discharge)",
    "A56.0":"Chlamydial infection of lower GU tract|Genitourinary|Chlamydia trachomatis|Doxycycline 100mg BD 7d or azithromycin 1g stat",
    "A56.1":"Chlamydial pelvic inflammatory disease|Genitourinary|Chlamydia trachomatis|Ceftriaxone + doxycycline + metronidazole",
    "A56.2":"Chlamydial infection of genitourinary tract, unspecified|Genitourinary|Chlamydia trachomatis|Doxycycline 100mg BD 7d",
    # ── MENINGITIS ──
    "G00.0":"Haemophilus meningitis|CNS|Haemophilus influenzae|Bacterial meningitis (community-acquired)",
    "G00.1":"Pneumococcal meningitis|CNS|Streptococcus pneumoniae|Bacterial meningitis (community-acquired)",
    "G00.2":"Streptococcal meningitis|CNS|Group B Streptococcus|Bacterial meningitis (community-acquired)",
    "G00.3":"Staphylococcal meningitis|CNS|Staphylococcus aureus|Vancomycin + rifampicin (post-neurosurgical)",
    "G00.8":"Other bacterial meningitis|CNS|Various|Bacterial meningitis (community-acquired)",
    "G00.9":"Bacterial meningitis, unspecified|CNS|Unknown|Bacterial meningitis (community-acquired)",
    "G01":"Meningitis in bacterial diseases classified elsewhere|CNS|TB, Listeria, syphilis|Investigate underlying cause",
    # ── PNEUMONIA ──
    "J13":"Pneumonia due to S. pneumoniae|Respiratory|Streptococcus pneumoniae|Community-acquired pneumonia (CAP)",
    "J14":"Pneumonia due to H. influenzae|Respiratory|Haemophilus influenzae|Community-acquired pneumonia (CAP)",
    "J15.0":"Pneumonia due to K. pneumoniae|Respiratory|Klebsiella pneumoniae|Hospital-acquired / Ventilator-associated pneumonia (HAP/VAP)",
    "J15.1":"Pneumonia due to Pseudomonas|Respiratory|Pseudomonas aeruginosa|Hospital-acquired / Ventilator-associated pneumonia (HAP/VAP)",
    "J15.2":"Pneumonia due to Staphylococcus|Respiratory|Staphylococcus aureus|Consider MRSA cover if risk factors",
    "J15.4":"Pneumonia due to other Streptococci|Respiratory|Streptococcus spp.|Community-acquired pneumonia (CAP)",
    "J15.5":"Pneumonia due to E. coli|Respiratory|Escherichia coli|Hospital-acquired / Ventilator-associated pneumonia (HAP/VAP)",
    "J15.6":"Pneumonia due to other Gram-negative bacteria|Respiratory|Gram-negative rods|Hospital-acquired / Ventilator-associated pneumonia (HAP/VAP)",
    "J15.9":"Bacterial pneumonia, unspecified|Respiratory|Unknown|Community-acquired pneumonia (CAP)",
    "J18.0":"Bronchopneumonia, unspecified|Respiratory|Mixed|Community-acquired pneumonia (CAP)",
    "J18.1":"Lobar pneumonia, unspecified|Respiratory|Streptococcus pneumoniae likely|Community-acquired pneumonia (CAP)",
    "J18.9":"Pneumonia, unspecified|Respiratory|Unknown|Community-acquired pneumonia (CAP)",
    "J85.1":"Lung abscess with pneumonia|Respiratory|Anaerobes, S. aureus, Klebsiella|Amoxicillin-clavulanate or clindamycin + drainage",
    "J86.9":"Empyema|Respiratory|S. pneumoniae, S. aureus, anaerobes|Chest drain + ceftriaxone + metronidazole",
    # ── ENT ──
    "H66.0":"Acute suppurative otitis media|ENT|S. pneumoniae, H. influenzae, Moraxella|Amoxicillin first-line; amoxicillin-clavulanate if fails",
    "H66.9":"Otitis media, unspecified|ENT|S. pneumoniae, H. influenzae|Amoxicillin PO; watchful waiting if mild",
    "J01.9":"Acute sinusitis, unspecified|ENT|S. pneumoniae, H. influenzae|Amoxicillin-clavulanate if bacterial suspected (>10d symptoms)",
    "J02.0":"Streptococcal pharyngitis|ENT|Group A Streptococcus|Penicillin V or amoxicillin; prevents rheumatic fever",
    "J03.0":"Streptococcal tonsillitis|ENT|Group A Streptococcus|Penicillin V or amoxicillin",
    "J36":"Peritonsillar abscess (quinsy)|ENT|Group A Strep, anaerobes, S. aureus|I&D + amoxicillin-clavulanate or clindamycin",
    # ── ENDOCARDITIS ──
    "I33.0":"Acute infective endocarditis|Cardiovascular|S. aureus, Streptococci, Enterococcus|Infective endocarditis",
    "I33.9":"Acute endocarditis, unspecified|Cardiovascular|S. aureus, Streptococci|Infective endocarditis",
    # ── URINARY TRACT ──
    "N10":"Acute pyelonephritis|Genitourinary|Escherichia coli, Klebsiella|Complicated UTI / Pyelonephritis",
    "N11.0":"Nonobstructive chronic pyelonephritis|Genitourinary|E. coli, Proteus, Klebsiella|Complicated UTI / Pyelonephritis",
    "N30.0":"Acute cystitis|Genitourinary|Escherichia coli|Uncomplicated UTI (cystitis)",
    "N30.9":"Cystitis, unspecified|Genitourinary|E. coli, Klebsiella, Enterococcus|Uncomplicated UTI (cystitis)",
    "N34.1":"Nonspecific urethritis|Genitourinary|Chlamydia, Mycoplasma, Ureaplasma|Doxycycline 100mg BD 7d",
    "N39.0":"UTI, site not specified|Genitourinary|Escherichia coli|Uncomplicated UTI (cystitis)",
    "N41.0":"Acute prostatitis|Genitourinary|E. coli, Klebsiella, Pseudomonas|Ciprofloxacin or TMP-SMX 4-6 weeks",
    # ── SKIN AND SOFT TISSUE ──
    "L01.0":"Impetigo|Skin|S. aureus, Group A Strep|Topical mupirocin; oral flucloxacillin if extensive",
    "L02.0":"Cutaneous abscess|Skin|Staphylococcus aureus|Skin and soft tissue infection (SSTI)",
    "L02.2":"Furuncle|Skin|Staphylococcus aureus|I&D primary; antibiotics if cellulitis/systemic signs",
    "L02.9":"Cutaneous abscess, unspecified|Skin|Staphylococcus aureus|Skin and soft tissue infection (SSTI)",
    "L03.0":"Cellulitis of finger/toe|Skin|S. aureus, Group A Strep|Flucloxacillin or amoxicillin-clavulanate",
    "L03.1":"Cellulitis of limb|Skin|Group A Strep, S. aureus|Skin and soft tissue infection (SSTI)",
    "L03.3":"Cellulitis of trunk|Skin|Group A Strep, S. aureus|Skin and soft tissue infection (SSTI)",
    "L03.9":"Cellulitis, unspecified|Skin|Group A Strep, S. aureus|Skin and soft tissue infection (SSTI)",
    "L08.0":"Pyoderma|Skin|S. aureus, Streptococcus|Flucloxacillin or amoxicillin-clavulanate",
    "L08.9":"Local infection of skin, unspecified|Skin|S. aureus|Skin and soft tissue infection (SSTI)",
    # ── BONE AND JOINT ──
    "M00.0":"Staphylococcal arthritis|Musculoskeletal|Staphylococcus aureus|Flucloxacillin IV (or vancomycin if MRSA); joint washout",
    "M00.9":"Pyogenic arthritis, unspecified|Musculoskeletal|S. aureus, Streptococci, GNR|Flucloxacillin + ceftriaxone pending culture; washout",
    "M86.0":"Acute haematogenous osteomyelitis|Musculoskeletal|S. aureus, Salmonella (sickle cell)|Flucloxacillin IV 4-6 weeks; clindamycin if MRSA",
    "M86.1":"Other acute osteomyelitis|Musculoskeletal|S. aureus, GNR|Culture-directed; surgical debridement",
    "M86.9":"Osteomyelitis, unspecified|Musculoskeletal|S. aureus|Flucloxacillin IV; prolonged course",
    # ── ABDOMINAL / SURGICAL ──
    "K35.2":"Acute appendicitis with peritonitis|Abdominal|E. coli, Bacteroides, Enterococcus|Intra-abdominal infection",
    "K35.9":"Acute appendicitis, unspecified|Abdominal|E. coli, Bacteroides|Intra-abdominal infection",
    "K57.2":"Diverticulitis with perforation/abscess|Abdominal|E. coli, Bacteroides fragilis|Intra-abdominal infection",
    "K61.0":"Anal abscess|Abdominal|E. coli, Bacteroides, Streptococci|I&D + metronidazole + ciprofloxacin if systemic",
    "K65.0":"Generalized acute peritonitis|Abdominal|E. coli, Bacteroides, Klebsiella|Intra-abdominal infection",
    "K65.9":"Peritonitis, unspecified|Abdominal|E. coli, anaerobes|Intra-abdominal infection",
    "K75.0":"Liver abscess|Abdominal|E. coli, Klebsiella, anaerobes, Entamoeba|Drainage + ceftriaxone + metronidazole",
    "K80.0":"Cholecystitis with cholelithiasis|Abdominal|E. coli, Klebsiella, Enterococcus|Ceftriaxone + metronidazole; cholecystectomy",
    "K83.0":"Cholangitis|Abdominal|E. coli, Klebsiella, Enterococcus|Piperacillin-tazobactam; ERCP if obstructed",
    # ── NEONATAL ──
    "P36.0":"Sepsis of newborn due to Group B Strep|Neonatal|Group B Streptococcus|Neonatal sepsis (early-onset)",
    "P36.1":"Sepsis of newborn due to other Streptococci|Neonatal|Streptococcus spp.|Neonatal sepsis (early-onset)",
    "P36.2":"Sepsis of newborn due to S. aureus|Neonatal|Staphylococcus aureus|Flucloxacillin + gentamicin; vancomycin if MRSA",
    "P36.3":"Sepsis of newborn due to other GNR|Neonatal|E. coli, Klebsiella|Neonatal sepsis (early-onset)",
    "P36.4":"Sepsis of newborn due to E. coli|Neonatal|Escherichia coli|Neonatal sepsis (early-onset)",
    "P36.9":"Bacterial sepsis of newborn, unspecified|Neonatal|GBS, E. coli, Klebsiella|Neonatal sepsis (early-onset)",
    "P23.0":"Congenital pneumonia due to viral agent|Neonatal|Viral|Supportive; add antibiotics if bacterial superinfection",
    "P23.4":"Congenital pneumonia due to E. coli|Neonatal|Escherichia coli|Ampicillin + gentamicin",
    # ── OBSTETRIC ──
    "O85":"Puerperal sepsis|Obstetric|Group A Strep, E. coli, anaerobes|Ampicillin + gentamicin + metronidazole; EMERGENCY",
    "O86.0":"Infection of obstetric surgical wound|Obstetric|S. aureus, E. coli, anaerobes|Amoxicillin-clavulanate; wound care",
    # ── SURGICAL / TRAUMA ──
    "T79.3":"Post-traumatic wound infection|Surgical|S. aureus, Streptococcus, GNR|Skin and soft tissue infection (SSTI)",
    "T81.4":"Infection following a procedure|Surgical|S. aureus, coagulase-neg Staph, GNR|Culture-directed; consider device removal if implant",
    "T84.5":"Infection of internal prosthetic joint|Musculoskeletal|S. aureus, CoNS, Streptococci|Vancomycin + rifampicin; surgical revision",
    # ── CRYPTOCOCCAL ──
    "B45.1":"Cerebral cryptococcosis|CNS|Cryptococcus neoformans|Cryptococcal meningitis (HIV)",
    "B45.9":"Cryptococcosis, unspecified|Systemic|Cryptococcus neoformans|Cryptococcal meningitis (HIV)",
    # ── PID ──
    "N70.0":"Acute salpingitis|Genitourinary|N. gonorrhoeae, Chlamydia, anaerobes|Ceftriaxone + doxycycline + metronidazole",
    "N73.0":"Acute parametritis/pelvic cellulitis|Genitourinary|Mixed aerobic/anaerobic|Ceftriaxone + doxycycline + metronidazole",
}

# Parse the pipe-delimited ICD entries into structured records
ICD_PARSED = {}
for code, val in ICD_10.items():
    parts = val.split("|")
    ICD_PARSED[code] = {
        "description": parts[0] if len(parts) > 0 else "",
        "category": parts[1] if len(parts) > 1 else "",
        "pathogens": parts[2] if len(parts) > 2 else "",
        "linked_syndrome": parts[3] if len(parts) > 3 else "",
    }

ICD_CATEGORIES = sorted(set(v["category"] for v in ICD_PARSED.values()))


AMR_RISK_FACTORS = {
    "Prior antibiotic use (last 90 days)":       {"weight":3,"affects":"ESBL, MRSA, CRE, FQ-R, CDI","reference":"WHO GLASS 2025; Tacconelli et al. Lancet ID 2018"},
    "Hospitalisation (last 90 days)":            {"weight":3,"affects":"MRSA, ESBL, CRE, VRE, Acinetobacter","reference":"IDSA HAP/VAP guidelines 2016"},
    "ICU admission (current or recent)":         {"weight":3,"affects":"MDR Gram-negatives, MRSA, Acinetobacter","reference":"Vincent et al. JAMA 2009 (EPIC II)"},
    "Indwelling device (catheter, CVC, ETT)":    {"weight":2,"affects":"MRSA, CRE, Pseudomonas, Candida","reference":"Magill et al. NEJM 2014"},
    "Immunosuppression (HIV, transplant, chemo)":{"weight":2,"affects":"Opportunistic infections, MDR-TB, iNTS, Cryptococcus","reference":"WHO OI guidelines 2021"},
    "HIV infection (CD4 <200)":                  {"weight":2,"affects":"TB (including MDR), iNTS, Cryptococcus, PJP, toxoplasma","reference":"UNAIDS/WHO"},
    "Diabetes mellitus":                         {"weight":1,"affects":"MRSA, diabetic foot polymicrobial, UTI","reference":"Pearson-Stuttard et al. Lancet Diabetes 2022"},
    "Chronic kidney disease":                    {"weight":1,"affects":"Dose adjustment needed, aminoglycoside/vancomycin toxicity risk","reference":"KDIGO guidelines"},
    "Travel to high-resistance region":          {"weight":2,"affects":"ESBL colonisation, CRE, MDR typhoid","reference":"Arcilla et al. Lancet ID 2017"},
    "Residence in long-term care facility":      {"weight":2,"affects":"MRSA, VRE, ESBL, CDI","reference":"Pop-Vicas et al. CID 2020"},
    "Age <1 year or >65 years":                  {"weight":1,"affects":"Higher susceptibility to invasive disease, neonatal ESBL","reference":"GLASS 2025"},
    "Known MDRO colonisation/prior culture":     {"weight":3,"affects":"Specific organism identified — use prior AST to guide therapy","reference":"IDSA stewardship guidelines 2016"},
    "Prior MRSA infection/colonisation":         {"weight":3,"affects":"MRSA SSTI, bacteraemia","reference":"Liu et al. CID 2011 (IDSA MRSA guidelines)"},
    "Livestock/agricultural antibiotic exposure": {"weight":1,"affects":"ESBL, colistin-R (mcr-1), FQ-R via One Health pathway","reference":"Van Boeckel et al. Science 2019"},
}

# ── 1F. REFERENCES ──────────────────────────────────────────────────────
REFERENCES = [
    ("WHO GLASS Report 2025","Global Antimicrobial Resistance and Use Surveillance System. WHO, 2025.","https://www.who.int/publications/i/item/9789240116337"),
    ("WHO AWaRe Classification 2023","Moja L, Zanichelli V, Mertz D, et al. Clin Microbiol Infect. 2024;30(S2):S1-S51.","https://doi.org/10.1016/j.cmi.2024.02.003"),
    ("MRSA in Africa — meta-analysis 2025","Proportion and antibiogram of MRSA in Africa. BMC Infect Dis. 2025.","https://doi.org/10.1186/s12879-025-10819-4"),
    ("East Africa AMR review","Tadesse BT, et al. A review of AMR in East Africa. Afr J Lab Med. 2017;6(1):346.","https://doi.org/10.4102/ajlm.v6i1.346"),
    ("AMR genetic diversity E. Africa","Salgado-Gaitán R, et al. Genetic diversity and risk factors for AMR transmission in East Africa. Antimicrob Resist Infect Control. 2020;9:127.","https://doi.org/10.1186/s13756-020-00786-3"),
    ("AMR mortality East Africa 2025","Adebisi YA, et al. Mortality burden of bacterial AMR in East Africa. Trop Med Health. 2025;53:19.","https://doi.org/10.1186/s41182-025-00870-x"),
    ("Surviving Sepsis 2021","Evans L, et al. Surviving Sepsis Campaign guidelines 2021. Crit Care Med. 2021;49(11):e1063-e1143.","https://doi.org/10.1097/CCM.0000000000005337"),
    ("IDSA MRSA Guidelines","Liu C, et al. MRSA infections in adults and children. CID. 2011;52(3):e18-e55.","https://doi.org/10.1093/cid/ciq146"),
    ("WHO Cryptococcal Meningitis 2022","Guidelines for diagnosing, preventing, and managing cryptococcal disease among adults, adolescents and children living with HIV. WHO, 2022.","https://www.who.int/publications/i/item/9789240052178"),
    ("Uganda Clinical Guidelines 2023","Ministry of Health Uganda. Uganda Clinical Guidelines 2023: National Guidelines for Management of Common Health Conditions. Kampala: MoH; 2023.","https://health.go.ug/"),
    ("Uganda AMR surveillance 2024","Muwanguzi D et al. Antibiotic Resistance in Select Tertiary Hospitals in Uganda 2020-2023. Antibiotics 2024;13(4):314.","https://doi.org/10.3390/antibiotics13040314"),
]

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                         2. LOGIC LAYER                                 ║
# ╚══════════════════════════════════════════════════════════════════════════╝

def get_resistance_level_color(pct):
    if pct >= 90: return "#2ecc71"
    elif pct >= 70: return "#f1c40f"
    elif pct >= 50: return "#e67e22"
    else: return "#e74c3c"

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║                         3. STREAMLIT UI                                ║
# ╚══════════════════════════════════════════════════════════════════════════╝

st.set_page_config(page_title="Rx Steward", page_icon="💊", layout="wide",
                   initial_sidebar_state="expanded")

st.markdown("""
<style>
    .stTabs [data-baseweb="tab-list"] { gap: 8px; }
    .stTabs [data-baseweb="tab"] { padding: 8px 16px; font-weight: 600; }
    div[data-testid="stMetricValue"] { font-size: 1.1rem; }
    .source-tag { font-size: 0.7em; color: #888; font-style: italic; }
</style>
""", unsafe_allow_html=True)

st.title("💊 Rx Steward")
st.caption("Antimicrobial Resistance & Stewardship Decision-Support Tool")

tabs = st.tabs(["🏥 Empiric Prescribing","🧫 Antibiogram Explorer",
                "⚠️ Risk Assessment","📖 Drug Library","📚 References"])

# ── TAB 1: EMPIRIC PRESCRIBING ──────────────────────────────────────────
# Risk factor → therapy modification rules
RISK_MODS = {
    "Prior MRSA infection/colonisation":         {"add":"+ Vancomycin IV","flag":"MRSA cover added based on colonisation/infection history."},
    "Known MDRO colonisation/prior culture":     {"add":"Escalate to broader-spectrum","flag":"Use prior AST to guide. Escalate empiric cover pending cultures."},
    "Prior antibiotic use (last 90 days)":       {"add":"Consider ESBL cover (carbapenem if serious)","flag":"Prior antibiotics increase ESBL, MRSA, and CDI risk. Avoid repeating same class."},
    "Hospitalisation (last 90 days)":            {"add":"Treat as healthcare-associated","flag":"Hospital pathogens (ESBL, MRSA, VRE) more likely. Use hospital-acquired regimens."},
    "ICU admission (current or recent)":         {"add":"+ Aminoglycoside (double Gram-neg cover); consider anti-pseudomonal","flag":"ICU pathogens: MDR Pseudomonas, Acinetobacter, ESBL. Double Gram-negative cover recommended."},
    "Indwelling device (catheter, CVC, ETT)":    {"add":"Consider device removal/change; add biofilm cover","flag":"Biofilm organisms: CoNS, Pseudomonas, Candida. Remove/change device before or with antibiotics."},
    "Immunosuppression (HIV, transplant, chemo)":{"add":"Broaden cover; consider opportunistic pathogens","flag":"Broader differential: PJP, TB, fungal, Nocardia. Add co-trimoxazole prophylaxis if not already on it."},
    "HIV infection (CD4 <200)":                  {"add":"TB workup; PJP cover; CrAg screen","flag":"High risk for TB, PJP, Cryptococcus, iNTS. Ensure co-trimoxazole prophylaxis. Screen CrAg."},
    "Diabetes mellitus":                         {"add":"Broaden to cover polymicrobial/MRSA","flag":"Diabetic patients: higher MRSA risk, polymicrobial foot infections. Ensure glycaemic control."},
    "Chronic kidney disease":                    {"add":"Dose-adjust renally cleared drugs","flag":"Avoid/adjust aminoglycosides, vancomycin (TDM essential). No nitrofurantoin if eGFR<30."},
    "Travel to high-resistance region":          {"add":"Consider ESBL colonisation","flag":"Post-travel ESBL colonisation rates 20-70%. Use carbapenem if seriously ill."},
    "Residence in long-term care facility":      {"add":"Treat as healthcare-associated","flag":"LTCF residents: higher MRSA, VRE, ESBL carriage. Use hospital-acquired pathogen spectrum."},
    "Age <1 year or >65 years":                  {"add":"Consider Listeria cover (ampicillin) if meningitis/sepsis","flag":"Age extremes: broader pathogen range. Add ampicillin for Listeria if CNS/sepsis in >50yr or neonates."},
    "Livestock/agricultural antibiotic exposure": {"add":"Consider ESBL/colistin-resistant organisms","flag":"One Health AMR pathway: ESBL, mcr-1 colistin resistance. Relevant for rural E. Africa."},
}

with tabs[0]:
    st.header("Empiric Prescribing Guide")

    # ── Build unified dropdown: syndromes + ICD-10 codes ──
    entry_options = []
    entry_map = {}
    # Add syndromes first
    for s_name in sorted(SYNDROMES.keys()):
        label = f"🏥 {s_name}"
        entry_options.append(label)
        entry_map[label] = {"type":"syndrome","syndrome":s_name}
    # Add ICD-10 codes grouped
    for code in sorted(ICD_PARSED.keys()):
        info = ICD_PARSED[code]
        linked = info["linked_syndrome"]
        if linked in SYNDROMES:
            label = f"📋 {code} — {info['description']}"
        else:
            label = f"📋 {code} — {info['description']} (guidance only)"
        entry_options.append(label)
        entry_map[label] = {"type":"icd","code":code,"info":info}

    selected_entry = st.selectbox("Search by clinical syndrome or ICD-10 code", entry_options,
                                   help="Type to search. Syndromes (🏥) have full protocols. ICD codes (📋) link to protocols where available.")

    region = st.selectbox("Region", ["Uganda","East Africa","Africa","Global"], key="emp_region")

    # Resolve to syndrome
    entry_data = entry_map[selected_entry]
    syn = None
    syndrome = None
    icd_code = None
    if entry_data["type"] == "syndrome":
        syndrome = entry_data["syndrome"]
        syn = SYNDROMES[syndrome]
    else:
        icd_code = entry_data["code"]
        icd_info = entry_data["info"]
        if icd_info["linked_syndrome"] in SYNDROMES:
            syndrome = icd_info["linked_syndrome"]
            syn = SYNDROMES[syndrome]
        else:
            st.subheader(f"{icd_code} — {icd_info['description']}")
            st.markdown(f"**Body system:** {icd_info['category']}  •  **Likely pathogens:** {icd_info['pathogens']}")
            st.info(f"**Empiric guidance:** {icd_info['linked_syndrome']}")
            st.caption("Full syndrome protocol not yet mapped for this code.")

    # ── Helper: get resistance justification ──
    def get_resistance_note(drug_name, pathogens, region):
        """Look up susceptibility of drug vs likely pathogens at given region."""
        notes = []
        for pathogen in pathogens:
            if pathogen in BACTERIA:
                bact_data = BACTERIA[pathogen]
                # Try region hierarchy
                for lvl in [region, "East Africa", "Africa", "Global"]:
                    if lvl in bact_data["susceptibility"]:
                        susc = bact_data["susceptibility"][lvl]
                        # Match drug name (handle combo drugs)
                        for abx_key, pct in susc.items():
                            if abx_key.lower() in drug_name.lower() or drug_name.lower().split()[0] in abx_key.lower():
                                emoji = "✅" if pct >= 85 else "⚠️" if pct >= 60 else "🚫"
                                short_name = pathogen.split("(")[0].strip()
                                if len(short_name) > 20:
                                    parts = short_name.split()
                                    short_name = parts[0][0] + ". " + " ".join(parts[1:])
                                notes.append(f"{emoji} {short_name}: {pct}% susceptible ({lvl})")
                                break
                        break
        return notes

    # ── Display therapy with justification ──
    if syn:
        st.subheader(syndrome)
        if icd_code:
            st.caption(f"Linked from ICD-10: {icd_code}")
        st.markdown(f"**Setting:** {syn['setting']}")
        st.markdown(f"**Likely pathogens:** {', '.join(syn['likely_pathogens'])}")

        aware_map = {"A":"🟢 Access","W":"🟡 Watch","R":"🔴 Reserve","A+W":"🟡 Access+Watch",
                     "W+W":"🟡 Watch","W+A":"🟡 Watch+Access","A+A":"🟢 Access","A+A+A":"🟢 Access",
                     "R/W":"🔴 Reserve/Watch","—":"⚪"}

        # Risk factors (visible, two columns)
        st.markdown("---")
        st.markdown("**Patient risk factors:**")
        selected_risks = []
        rf_cols = st.columns(2)
        rf_list = list(AMR_RISK_FACTORS.keys())
        for i, rf in enumerate(rf_list):
            with rf_cols[i % 2]:
                if st.checkbox(rf, key=f"emp2_rf_{rf}"):
                    selected_risks.append(rf)

        # Risk modifications
        mods = []
        for rf in selected_risks:
            if rf in RISK_MODS:
                mods.append(RISK_MODS[rf])

        st.markdown("---")

        # First-line with resistance justification
        if mods:
            st.markdown("### ⚡ Risk-adjusted empiric therapy")
            st.error("**Risk factors detected** — therapy modified. Obtain cultures BEFORE antibiotics.")
        else:
            st.markdown("### First-line empiric therapy")

        for drug, dose, aware in syn["first_line"]:
            aware_label = aware_map.get(aware, aware)
            st.markdown(f"**{drug}** — {dose}  `{aware_label}`")
            # Resistance justification
            res_notes = get_resistance_note(drug, syn["likely_pathogens"], region)
            if res_notes:
                st.caption("  ".join(res_notes))

        if mods:
            st.markdown("**Risk-adjusted additions:**")
            for m in mods:
                st.markdown(f"🔺 **{m['add']}** — {m['flag']}")

        # Second-line (visible, not collapsed)
        st.markdown("### Second-line / alternatives")
        for drug, dose, aware in syn["second_line"]:
            aware_label = aware_map.get(aware, aware)
            st.markdown(f"**{drug}** — {dose}  `{aware_label}`")
            res_notes = get_resistance_note(drug, syn["likely_pathogens"], region)
            if res_notes:
                st.caption("  ".join(res_notes))

        if syn.get("avoid"):
            st.warning(f"**Avoid / caution:** {syn['avoid']}")
        st.info(f"**Clinical notes:** {syn['notes']}")

# ── TAB 2: ANTIBIOGRAM EXPLORER ─────────────────────────────────────────
with tabs[1]:
    st.header("Antibiogram Explorer")
    col1, col2 = st.columns([1, 3])
    with col1:
        org = st.selectbox("Organism", list(BACTERIA.keys()))
        bact = BACTERIA[org]
        # Check for regional/facility data in session state
        all_levels = ["Global","Africa","East Africa","Uganda"]
        if "regional_data" in st.session_state and isinstance(st.session_state.get("regional_data"), dict):
            for rname, rdata in st.session_state["regional_data"].items():
                if isinstance(rdata, dict) and org in rdata and rname not in all_levels:
                    all_levels.append(rname)
        available_regions = [l for l in all_levels if l in bact["susceptibility"]]
        # Also check regional uploads
        if "regional_data" in st.session_state and isinstance(st.session_state.get("regional_data"), dict):
            for l in all_levels:
                if l not in available_regions and l in st.session_state["regional_data"] and org in st.session_state["regional_data"].get(l, {}):
                    available_regions.append(l)
        if not available_regions:
            available_regions = ["Global"]
        reg = st.selectbox("Region (most local shown first)", list(reversed(available_regions)), key="abg_region")
        avail = [l for l in all_levels if l in bact["susceptibility"]]
        st.caption(f"✅ Data: {', '.join(avail)}")
        st.markdown("---")
        with st.expander("📤 Upload facility antibiogram"):
            st.markdown("Upload a CSV with columns: `organism`, `antibiotic`, `percent_susceptible`")
            uploaded = st.file_uploader("Facility antibiogram CSV", type="csv", key="facility_csv")
            if uploaded:
                try:
                    fac_df = pd.read_csv(uploaded)
                    fac_dict = {}
                    for _, row in fac_df.iterrows():
                        o = str(row.get("organism","")).strip()
                        a = str(row.get("antibiotic","")).strip()
                        v = float(row.get("percent_susceptible",0))
                        if o not in fac_dict: fac_dict[o] = {}
                        fac_dict[o][a] = v
                    st.session_state["facility_data"] = fac_dict
                    st.success(f"Loaded {len(fac_dict)} organisms from facility antibiogram.")
                except Exception as e:
                    st.error(f"Error reading CSV: {e}")
    with col2:
        st.subheader(org)
        st.markdown(f"**Classification:** {bact['gram']} | {bact['family']}")
        st.markdown(f"**Key diseases:** {', '.join(bact['diseases'])}")
        st.markdown(f"**Resistance mechanisms:** {', '.join(bact['resistance_mechanisms'])}")
        st.info(bact["key_notes"])
        if "sources" in bact:
            with st.expander("📖 Data sources for this organism"):
                for level, src in bact["sources"].items():
                    st.markdown(f"**{level}:** {src}")

        if "regional_data" in st.session_state and isinstance(st.session_state.get("regional_data"), dict) and reg in st.session_state["regional_data"] and org in st.session_state["regional_data"][reg]:
            susc = st.session_state["regional_data"][reg][org]
        else:
            susc = bact["susceptibility"].get(reg, bact["susceptibility"].get("Global",{}))
        if susc:
            st.markdown(f"### Susceptibility profile — {reg}")
            df = pd.DataFrame([{"Antibiotic": k, "% Susceptible": v,
                                "AWaRe": ANTIBIOTICS.get(k,{}).get("aware","?"),
                                "Interpretation": "Reliable" if v>=85 else "Caution" if v>=60 else "Avoid empirically"}
                               for k,v in susc.items()])
            df = df.sort_values("% Susceptible", ascending=False)

            fig = px.bar(df, x="% Susceptible", y="Antibiotic", orientation="h",
                         color="% Susceptible", range_x=[0,105],
                         color_continuous_scale=["#e74c3c","#e67e22","#f1c40f","#2ecc71"],
                         range_color=[0,100], height=max(300, len(df)*30))
            fig.update_layout(yaxis=dict(autorange="reversed"), showlegend=False, margin=dict(l=0,r=0,t=20,b=0))
            fig.add_vline(x=85, line_dash="dash", line_color="green", annotation_text="≥85% = reliable empiric")
            fig.add_vline(x=60, line_dash="dot", line_color="orange", annotation_text="<60% = avoid empiric")
            st.plotly_chart(fig, width='stretch')
            st.dataframe(df, hide_index=True, width='stretch')

# ── TAB 4: RISK ASSESSMENT ──────────────────────────────────────────────
with tabs[2]:
    st.header("AMR Risk Assessment")
    st.markdown("Systematic evaluation of patient-level resistance risk factors to guide empiric therapy escalation and diagnostic prioritisation.")

    st.subheader("Risk Factors")
    amr_sel = []
    for rf, info in AMR_RISK_FACTORS.items():
        if st.checkbox(rf, key=f"risk_amr_{rf}"):
            amr_sel.append((rf, info))

    if amr_sel:
        amr_total = sum(i["weight"] for _,i in amr_sel)
        amr_level = "🟢 Standard empiric" if amr_total <=3 else "🟡 Consider broader spectrum" if amr_total <=6 else "🔴 High MDR risk — escalate empiric cover"
        st.metric("AMR risk level", amr_level)

        organisms_affected = set()
        for _,i in amr_sel:
            for o in i["affects"].split(","):
                organisms_affected.add(o.strip())
        st.markdown(f"**Organisms/resistance patterns to cover:** {', '.join(sorted(organisms_affected))}")

        st.markdown("---")
        st.subheader("Recommendation")
        if amr_total > 6:
            st.error("**High MDR risk.** This patient requires: (1) cultures/sensitivity BEFORE empiric antibiotics, (2) broader-spectrum empiric cover, (3) ID specialist input if available, (4) mandatory de-escalation review at 48-72h.")
        elif amr_total > 3:
            st.warning("**Moderate risk.** Consider broader-spectrum empiric cover with planned de-escalation. Obtain cultures before antibiotics where possible.")
        else:
            st.success("**Standard risk.** Guideline-based empiric therapy appropriate. Monitor clinical response at 48-72h.")

        st.markdown("### Selected risk factors")
        for rf, info in amr_sel:
            st.markdown(f"- **{rf}** — affects: {info['affects']}  \n  *Ref: {info['reference']}*")

# ── TAB 5: DRUG LIBRARY ─────────────────────────────────────────────────
with tabs[3]:
    st.header("Antibiotic Library")
    st.markdown(f"**{len(ANTIBIOTICS)} antibiotics** classified by WHO AWaRe (Access / Watch / Reserve).")

    aware_filter = st.multiselect("Filter by AWaRe", ["A","W","R"], default=["A","W","R"],
                                  format_func=lambda x:{"A":"🟢 Access","W":"🟡 Watch","R":"🔴 Reserve"}[x])
    class_filter = st.multiselect("Filter by class", sorted(set(a["class"] for a in ANTIBIOTICS.values())))

    for name, info in sorted(ANTIBIOTICS.items()):
        if info["aware"] not in aware_filter:
            continue
        if class_filter and info["class"] not in class_filter:
            continue
        aware_label = {"A":"🟢 Access","W":"🟡 Watch","R":"🔴 Reserve"}.get(info["aware"],"")
        with st.expander(f"{aware_label}  **{name}** — {info['class']}"):
            st.markdown(f"**Route:** {info['route']}  \n**Spectrum:** {info['spectrum']}  \n**Notes:** {info['notes']}")

# ── TAB 6: REFERENCES ───────────────────────────────────────────────────
with tabs[4]:
    st.header("Data Sources & References")
    st.markdown("All clinical data in this tool is derived from the following peer-reviewed and institutional sources.")
    for i, (title, citation, url) in enumerate(REFERENCES, 1):
        st.markdown(f"**{i}. {title}**  \n{citation}  \n[{url}]({url})")
    st.markdown("---")
    st.caption("Rx Steward v2.1 — aligned with Uganda Clinical Guidelines 2023 — For clinical decision support only. Not a substitute for clinical judgement, local antibiogram data, or culture-directed therapy. Always de-escalate based on susceptibility results.")
