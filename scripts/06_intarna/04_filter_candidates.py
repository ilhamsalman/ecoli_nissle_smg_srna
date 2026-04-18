import os
import csv
import pandas as pd
from collections import defaultdict

# ─────────────────────────────────────────────
# Load data
# ─────────────────────────────────────────────
base_dir = os.path.join("results", "sRNA_detection", "06_intarna")

exp_df = pd.read_csv(os.path.join(base_dir, "exp_predictions_combined.csv"), sep=";")
sta_df = pd.read_csv(os.path.join(base_dir, "sta_predictions_combined.csv"), sep=";")

# Rename energy column
exp_df.rename(columns={"E": "deltaG"}, inplace=True)
sta_df.rename(columns={"E": "deltaG"}, inplace=True)

# Drop header-duplicate rows if any
exp_df = exp_df[exp_df["srna"] != "srna"]
sta_df = sta_df[sta_df["srna"] != "srna"]

# Convert deltaG to float
exp_df["deltaG"] = pd.to_numeric(exp_df["deltaG"], errors="coerce")
sta_df["deltaG"] = pd.to_numeric(sta_df["deltaG"], errors="coerce")

# Remove N_A targets (unnamed/unannotated)
exp_named = exp_df[~exp_df["target"].str.startswith("N_A")].copy()
sta_named = sta_df[~sta_df["target"].str.startswith("N_A")].copy()

print(f"EXP total interactions: {len(exp_df)}, named-target interactions: {len(exp_named)}")
print(f"STA total interactions: {len(sta_df)}, named-target interactions: {len(sta_named)}")

# ─────────────────────────────────────────────
# Gene functional annotation dictionary
# ─────────────────────────────────────────────
gene_functions = {
    # Stress response / chaperones
    "cspA": ("Cold shock protein A", "Stress response"),
    "htpG": ("High-temperature protein G (Hsp90 homolog)", "Heat stress / chaperone"),
    "hscA": ("Heat shock cognate protein HscA (Hsp70)", "Iron-sulfur cluster assembly / chaperone"),
    "hdeA": ("Acid resistance protein HdeA", "Acid stress response"),
    "hdeB": ("Acid resistance protein HdeB", "Acid stress response"),
    "hdeD": ("Acid resistance inner membrane protein", "Acid stress response"),
    "gadA": ("Glutamate decarboxylase A", "Acid stress response (GDAR system)"),
    "gadA_2": ("Glutamate decarboxylase A (isoform)", "Acid stress response (GDAR system)"),
    "gadB": ("Glutamate decarboxylase B", "Acid stress response (GDAR system)"),
    "gadB_2": ("Glutamate decarboxylase B (isoform)", "Acid stress response (GDAR system)"),
    "slp": ("Outer membrane lipoprotein Slp", "Stationary phase / stress response"),
    "osmC": ("Osmotically inducible protein OsmC", "Osmotic stress response"),
    "osmF": ("Osmotically inducible protein OsmF", "Osmotic stress response"),
    "ibpA": ("Small heat shock protein IbpA", "Heat/oxidative stress"),
    "ibpB": ("Small heat shock protein IbpB", "Heat/oxidative stress"),
    "hspQ": ("Heat shock protein HspQ", "Protein quality control / stress"),
    "clpP": ("ATP-dependent Clp protease proteolytic subunit", "Protein degradation / stress"),
    "clpX": ("ATP-dependent Clp protease ATPase subunit ClpX", "Protein degradation / stress"),
    "hslU": ("ATP-dependent protease ATPase subunit HslU", "Protein degradation / stress"),
    "htpX": ("Heat shock protein HtpX", "Membrane protein quality control"),
    "dsbG": ("Disulfide bond isomerase G", "Oxidative stress / periplasm"),
    "grxA": ("Glutaredoxin 1", "Oxidative stress response"),
    "trxC": ("Thioredoxin 2", "Oxidative stress response"),
    "ybfE": ("Predicted oxidoreductase", "Oxidative stress response"),
    "elaB": ("Stress-induced protein ElaB", "Membrane stress / stationary phase"),

    # Iron / metal homeostasis
    "feoB": ("Ferrous iron transporter B", "Iron acquisition"),
    "znuA": ("Zinc ABC transporter periplasmic binding protein", "Zinc homeostasis"),
    "znuB": ("Zinc ABC transporter membrane component", "Zinc homeostasis"),
    "znuC": ("Zinc ABC transporter ATP-binding component", "Zinc homeostasis"),
    "zinT": ("Zinc-binding periplasmic protein ZinT", "Zinc stress response"),
    "copA": ("Copper-transporting ATPase CopA", "Copper homeostasis"),
    "cueO": ("Multicopper oxidase CueO", "Copper homeostasis"),
    "cusA": ("Copper/silver efflux pump CusA", "Copper/silver resistance"),
    "cusB": ("Copper/silver efflux pump CusB", "Copper/silver resistance"),
    "rcnA": ("Nickel/cobalt efflux protein RcnA", "Nickel/cobalt resistance"),
    "iscR": ("Iron-sulfur cluster regulator IscR", "Iron-sulfur cluster assembly"),
    "iscA": ("Iron-sulfur cluster assembly protein IscA", "Iron-sulfur cluster assembly"),
    "iscU": ("Iron-sulfur cluster assembly scaffold IscU", "Iron-sulfur cluster assembly"),
    "zraP": ("Zinc resistance-associated protein ZraP", "Zinc stress response"),
    "modA": ("Molybdate-binding periplasmic protein ModA", "Molybdenum transport"),
    "modC": ("Molybdate ABC transporter ATP-binding protein", "Molybdenum transport"),
    "fhuA": ("Ferrichrome outer membrane transporter FhuA", "Iron acquisition (siderophore)"),
    "fiu": ("Ferric iron uptake protein Fiu", "Iron acquisition"),
    "fyuA": ("Yersiniabactin receptor FyuA", "Siderophore-mediated virulence"),
    "ybtP": ("Yersiniabactin ABC transporter", "Iron acquisition / virulence"),
    "ybtQ": ("Yersiniabactin ABC transporter", "Iron acquisition / virulence"),
    "ybtX": ("Yersiniabactin biosynthesis protein", "Iron acquisition / virulence"),
    "irp9": ("Yersiniabactin biosynthesis enzyme Irp9", "Iron acquisition / virulence"),
    "HMWP2": ("High-molecular-weight protein 2 (yersiniabactin)", "Iron acquisition / virulence"),

    # Motility / chemotaxis
    "fliC": ("Flagellin FliC", "Flagellar motility"),
    "cheZ": ("Chemotaxis phosphatase CheZ", "Chemotaxis / motility"),
    "ycgR": ("Flagellar brake protein YcgR", "c-di-GMP-mediated motility control"),

    # Biofilm / surface structures
    "fimA": ("Major type 1 fimbrial subunit FimA", "Biofilm / adhesion"),
    "fimC": ("Fimbrial chaperone FimC", "Biofilm / adhesion"),
    "fimD": ("Fimbrial usher protein FimD", "Biofilm / adhesion"),
    "fimG": ("Fimbrial protein FimG", "Biofilm / adhesion"),
    "fimI": ("Fimbrial protein FimI", "Biofilm / adhesion"),
    "ecpA": ("E. coli common pilus subunit EcpA", "Biofilm / adhesion"),
    "ecpB": ("E. coli common pilus chaperone EcpB", "Biofilm / adhesion"),
    "ecpC": ("E. coli common pilus usher EcpC", "Biofilm / adhesion"),
    "ecpR": ("E. coli common pilus regulator EcpR", "Biofilm / adhesion"),
    "csgA": ("Major curlin subunit CsgA", "Biofilm / curli"),
    "csgB": ("Minor curlin subunit CsgB", "Biofilm / curli"),
    "bhsA": ("Biofilm-related outer membrane protein BhsA", "Biofilm"),
    "bssS": ("Biofilm formation regulator BssS", "Biofilm"),
    "bluF": ("Blue light-sensing anti-repressor BluF", "Biofilm / c-di-GMP signaling"),

    # Metabolism
    "adhE": ("Alcohol/aldehyde dehydrogenase AdhE", "Anaerobic fermentation"),
    "aspA": ("Aspartate ammonia-lyase AspA", "Amino acid metabolism"),
    "fumB": ("Fumarate hydratase B (anaerobic)", "TCA cycle / anaerobic metabolism"),
    "fumB_2": ("Fumarate hydratase B (isoform)", "TCA cycle / anaerobic metabolism"),
    "gapA": ("Glyceraldehyde-3-phosphate dehydrogenase", "Glycolysis"),
    "pfkA": ("6-phosphofructokinase I", "Glycolysis"),
    "mqo": ("Malate:quinone oxidoreductase", "TCA cycle"),
    "sdhA": ("Succinate dehydrogenase flavoprotein subunit", "TCA cycle"),
    "cyoA": ("Cytochrome bo3 ubiquinol oxidase subunit I", "Aerobic respiration"),
    "nuoA": ("NADH:ubiquinone oxidoreductase subunit A", "Respiratory chain"),
    "nuoB": ("NADH:ubiquinone oxidoreductase subunit B", "Respiratory chain"),
    "nuoN": ("NADH:ubiquinone oxidoreductase subunit N", "Respiratory chain"),
    "dcuB": ("Anaerobic C4-dicarboxylate transporter DcuB", "Anaerobic metabolism"),
    "focA": ("Formate channel FocA", "Formate export / anaerobic metabolism"),
    "dadX": ("Alanine racemase DadX", "Amino acid metabolism"),
    "leuC": ("3-isopropylmalate dehydratase large subunit", "Leucine biosynthesis"),
    "lysA": ("Diaminopimelate decarboxylase LysA", "Lysine biosynthesis"),
    "betA": ("Choline dehydrogenase BetA", "Osmoprotectant synthesis"),
    "betB": ("Betaine aldehyde dehydrogenase BetB", "Osmoprotectant synthesis"),
    "betI": ("Bet operon transcriptional repressor BetI", "Osmoprotectant regulation"),
    "betT": ("Choline transporter BetT", "Osmoprotectant synthesis"),
    "tauA": ("Taurine ABC transporter periplasmic protein", "Sulfur acquisition"),
    "sbp": ("Sulfate-binding protein Sbp", "Sulfur acquisition"),
    "tcyP": ("Cystine transporter TcyP", "Sulfur/cysteine acquisition"),
    "tcyJ": ("Cystine ABC transporter component", "Sulfur/cysteine acquisition"),
    "tcyN": ("Cystine ABC transporter component", "Sulfur/cysteine acquisition"),
    "bioD": ("Dethiobiotin synthase BioD", "Biotin biosynthesis"),
    "cobW": ("Cobalamin biosynthesis protein CobW", "Cobalamin/B12 synthesis"),
    "gdhA": ("Glutamate dehydrogenase GdhA", "Nitrogen/amino acid metabolism"),
    "metA": ("Homoserine O-succinyltransferase MetA", "Methionine biosynthesis"),
    "metK": ("S-adenosylmethionine synthetase MetK", "Methionine metabolism / SAM"),
    "metQ": ("Methionine ABC transporter periplasmic protein", "Methionine transport"),
    "argB": ("N-acetylglutamate kinase ArgB", "Arginine biosynthesis"),
    "argO": ("Arginine exporter ArgO", "Arginine export"),
    "carA": ("Carbamoyl phosphate synthase small subunit", "Pyrimidine/arginine biosynthesis"),
    "carB": ("Carbamoyl phosphate synthase large subunit", "Pyrimidine/arginine biosynthesis"),
    "serC": ("Phosphoserine aminotransferase SerC", "Serine biosynthesis"),
    "thrC": ("Threonine synthase ThrC", "Threonine biosynthesis"),
    "guaA": ("GMP synthetase GuaA", "Purine biosynthesis"),
    "nadB": ("L-aspartate oxidase NadB", "NAD biosynthesis"),
    "nadE": ("NAD+ synthetase NadE", "NAD biosynthesis"),
    "pncB": ("Nicotinate phosphoribosyltransferase PncB", "NAD salvage"),
    "ribB": ("3,4-dihydroxy-2-butanone 4-phosphate synthase", "Riboflavin biosynthesis"),
    "menF": ("Isochorismate synthase MenF", "Menaquinone biosynthesis"),
    "ispG": ("4-hydroxy-3-methylbut-2-enyl diphosphate synthase", "Isoprenoid biosynthesis"),
    "tktB": ("Transketolase B", "Pentose phosphate pathway"),
    "edd": ("Phosphogluconate dehydratase Edd", "Entner-Doudoroff pathway"),
    "gntK": ("Gluconate kinase GntK", "Gluconate metabolism"),
    "manA": ("Mannose-6-phosphate isomerase ManA", "Sugar metabolism"),
    "pgl": ("6-phosphogluconolactonase Pgl", "Pentose phosphate pathway"),
    "prsA": ("Ribose-phosphate pyrophosphokinase PrsA", "Nucleotide biosynthesis"),
    "aroA": ("3-phosphoshikimate 1-carboxyvinyltransferase", "Aromatic amino acid biosynthesis"),
    "hisG": ("ATP phosphoribosyltransferase HisG", "Histidine biosynthesis"),
    "glmM": ("Phosphoglucosamine mutase GlmM", "Peptidoglycan precursor synthesis"),
    "rfaG": ("Lipopolysaccharide core biosynthesis protein RfaG", "LPS biosynthesis / envelope"),
    "rfaP": ("LPS core heptose kinase RfaP", "LPS biosynthesis / envelope"),
    "rfaQ": ("LPS core biosynthesis protein RfaQ", "LPS biosynthesis / envelope"),
    "lpxP": ("Palmitoleoyl-ACP-dependent acyltransferase LpxP", "LPS/lipid A biosynthesis"),
    "kfiA": ("Polysialic acid biosynthesis protein KfiA", "Capsule biosynthesis"),
    "kfiA_2": ("Polysialic acid biosynthesis protein KfiA (isoform)", "Capsule biosynthesis"),
    "kfiB": ("Polysialic acid biosynthesis protein KfiB", "Capsule biosynthesis"),
    "kfiB_2": ("Polysialic acid biosynthesis protein KfiB (isoform)", "Capsule biosynthesis"),
    "kfiB_3": ("Polysialic acid biosynthesis protein KfiB (isoform 2)", "Capsule biosynthesis"),

    # Virulence
    "ompT": ("Outer membrane protease OmpT", "Virulence / proteolysis"),
    "clbC": ("Colibactin biosynthesis protein ClbC", "Genotoxin / virulence"),
    "clbJ": ("Colibactin biosynthesis protein ClbJ", "Genotoxin / virulence"),
    "icsS": ("IcsS protein", "Virulence / actin motility"),

    # Envelope / outer membrane
    "yhbU": ("Predicted peptidase YhbU", "Envelope stress / peptidoglycan"),
    "yhbV": ("Predicted peptidase YhbV", "Envelope stress / peptidoglycan"),
    "pepE": ("Peptidase E", "Peptide degradation"),
    "pepT": ("Peptidase T", "Peptide degradation"),
    "nlpA": ("Lipoprotein NlpA", "Outer membrane / envelope"),
    "tatB": ("Twin-arginine translocase subunit TatB", "Protein secretion / envelope"),
    "tatC": ("Twin-arginine translocase subunit TatC", "Protein secretion / envelope"),
    "mepM": ("Murein endopeptidase MepM", "Peptidoglycan remodeling"),
    "lysM": ("Peptidoglycan-binding protein LysM", "Cell wall binding"),

    # DNA repair / replication
    "recF": ("DNA repair protein RecF", "DNA repair / SOS response"),
    "recQ": ("ATP-dependent DNA helicase RecQ", "DNA repair"),
    "uvrD": ("DNA helicase II UvrD", "DNA repair"),
    "priA": ("Primosome assembly protein PriA", "DNA replication restart"),
    "priA_2": ("Primosome assembly protein PriA (isoform)", "DNA replication restart"),

    # Transcription / translation / regulation
    "stpA": ("DNA-binding protein StpA (H-NS paralog)", "Nucleoid-associated protein / gene regulation"),
    "rplL": ("50S ribosomal protein L7/L12", "Translation"),
    "rpmE": ("50S ribosomal protein L31", "Translation / zinc-responsive"),
    "rpmF": ("50S ribosomal protein L32", "Translation"),
    "rpsT": ("30S ribosomal protein S20", "Translation"),
    "rpsU": ("30S ribosomal protein S21", "Translation"),
    "rsgA": ("Ribosome small subunit-dependent GTPase A", "Ribosome biogenesis"),
    "rimM": ("Ribosome maturation factor RimM", "Ribosome biogenesis"),
    "hpf": ("Hibernation-promoting factor Hpf", "Ribosome hibernation / stationary phase"),
    "dbpA": ("ATP-dependent RNA helicase DbpA", "Ribosome biogenesis"),
    "tufA": ("Elongation factor Tu (EF-Tu)", "Translation / stress"),
    "pcnB": ("Poly(A) polymerase PcnB", "mRNA degradation / regulation"),
    "rbsR": ("Ribose operon repressor RbsR", "Sugar transport regulation"),
    "aceK": ("Isocitrate dehydrogenase kinase/phosphatase AceK", "TCA cycle regulation"),
    "cbl": ("Sulfate regulon transcriptional activator Cbl", "Sulfur metabolism regulation"),
    "rstA": ("Response regulator RstA", "Two-component signal transduction"),
    "rstB": ("Sensor histidine kinase RstB", "Two-component signal transduction"),
    "alx": ("Predicted inner membrane protein Alx", "Alkaline stress response"),
    "btsT": ("Pyruvate sensor/transporter BtsT", "Carbon source sensing"),
    "caiF": ("Carnitine metabolism transcriptional activator CaiF", "Carnitine metabolism"),
    "caiC": ("Carnitine CoA-ligase CaiC", "Carnitine metabolism"),

    # Miscellaneous / SMG-relevant
    "glpE": ("Thiosulfate sulfurtransferase GlpE", "Sulfur metabolism"),
    "yfcC": ("Predicted inner membrane protein YfcC", "Unknown / membrane"),
    "ycgK": ("Predicted inner membrane protein YcgK", "Unknown / membrane"),
    "yacL": ("Predicted membrane protein YacL", "Unknown / membrane"),
    "frmA": ("S-(hydroxymethyl)glutathione dehydrogenase FrmA", "Formaldehyde detoxification"),
    "solA": ("N-methyl amino acid oxidase SolA", "Oxidative metabolism"),
    "fadE": ("Acyl-CoA dehydrogenase FadE", "Fatty acid beta-oxidation"),
    "deoD": ("Purine nucleoside phosphorylase DeoD", "Nucleoside catabolism"),
    "entH": ("Enterobactin biosynthesis thioesterase EntH", "Siderophore biosynthesis"),
    "gpmM": ("Phosphoglycerate mutase GpmM", "Glycolysis"),
    "coaE": ("Dephospho-CoA kinase CoaE", "Coenzyme A biosynthesis"),
    "lipB": ("Octanoyltransferase LipB", "Lipoic acid biosynthesis"),
    "cspD": ("Cold shock protein D", "Stress response / growth arrest"),
    "nfeR": ("Nitrite/nitrate response regulator NfeR", "Nitrogen metabolism"),
    "narZ": ("Nitrate reductase Z alpha subunit", "Anaerobic nitrate respiration"),
    "ghrB": ("Glyoxylate/hydroxypyruvate reductase B", "Glyoxylate metabolism"),
    "glcD": ("Glycolate oxidase subunit GlcD", "Glycolate metabolism"),
    "hycA": ("Hydrogenase-3 maturation factor HycA", "Hydrogen metabolism"),
    "hycB": ("Hydrogenase-3 iron-sulfur protein HycB", "Hydrogen metabolism"),
    "hycC": ("Hydrogenase-3 membrane subunit HycC", "Hydrogen metabolism"),
    "hycD": ("Hydrogenase-3 membrane subunit HycD", "Hydrogen metabolism"),
    "hycE": ("Hydrogenase-3 large subunit HycE", "Hydrogen metabolism"),
    "hycH": ("Hydrogenase-3 maturation protein HycH", "Hydrogen metabolism"),
    "hypB": ("Hydrogenase maturation factor HypB", "Hydrogenase maturation"),
    "hypF": ("Hydrogenase maturation factor HypF", "Hydrogenase maturation"),
    "fdhD": ("Formate dehydrogenase accessory protein FdhD", "Formate metabolism"),
    "adiA": ("Biodegradative arginine decarboxylase AdiA", "Acid resistance / arginine metabolism"),
    "gltX": ("Glutamyl-tRNA synthetase GltX", "Translation / aminoacyl-tRNA"),
    "gltX_2": ("Glutamyl-tRNA synthetase GltX (isoform)", "Translation / aminoacyl-tRNA"),
    "lysS": ("Lysyl-tRNA synthetase LysS", "Translation / aminoacyl-tRNA"),
    "ppa": ("Inorganic pyrophosphatase Ppa", "Energy metabolism"),
    "hcyG": ("Homocysteine S-methyltransferase HcyG", "Methionine cycle"),
    "prmB": ("Ribosomal protein L3 methyltransferase PrmB", "Ribosome modification"),
    "queC": ("7-cyano-7-deazaguanine synthase QueC", "tRNA modification"),
    "thiM": ("Hydroxyethylthiazole kinase ThiM", "Thiamine biosynthesis"),
    "tppB": ("Thiamine pyrophosphate transporter TppB", "Thiamine transport"),
    "dtpD": ("Dipeptide/tripeptide permease DtpD", "Peptide transport"),
    "yejM": ("Lipid A modification enzyme YejM", "Envelope / LPS modification"),
    "zapD": ("Cell division protein ZapD", "Cell division"),
    "kbp": ("Klebsiella-type bacteriocin precursor Kbp", "Bacteriocin / competition"),
}

# ─────────────────────────────────────────────
# Priority categories
# ─────────────────────────────────────────────
priority_categories = {
    "Stress response", "Heat stress / chaperone", "Acid stress response (GDAR system)",
    "Acid stress response", "Osmotic stress response", "Oxidative stress response",
    "Oxidative stress / periplasm", "Protein degradation / stress", "Membrane protein quality control",
    "Protein quality control / stress", "Heat/oxidative stress", "Membrane stress / stationary phase",
    "Iron acquisition", "Zinc homeostasis", "Zinc stress response", "Copper homeostasis",
    "Copper/silver resistance", "Iron-sulfur cluster assembly", "Iron acquisition (siderophore)",
    "Iron acquisition / virulence", "Siderophore-mediated virulence",
    "Flagellar motility", "Chemotaxis / motility", "c-di-GMP-mediated motility control",
    "Biofilm / adhesion", "Biofilm / curli", "Biofilm", "Biofilm / c-di-GMP signaling",
    "Virulence / proteolysis", "Genotoxin / virulence", "Virulence / actin motility",
    "Anaerobic fermentation", "TCA cycle / anaerobic metabolism", "Anaerobic metabolism",
    "Respiratory chain", "Aerobic respiration", "TCA cycle",
    "Osmoprotectant synthesis", "Osmoprotectant regulation",
    "LPS biosynthesis / envelope", "Capsule biosynthesis", "Envelope stress / peptidoglycan",
    "Protein secretion / envelope", "Peptidoglycan remodeling",
    "DNA repair / SOS response", "DNA repair",
    "Ribosome hibernation / stationary phase",
    "Nucleoid-associated protein / gene regulation",
    "Two-component signal transduction",
    "Alkaline stress response", "Cold shock / stress",
    "Hydrogen metabolism", "Formate metabolism", "Anaerobic nitrate respiration",
    "Fatty acid beta-oxidation", "Siderophore biosynthesis",
}

def annotate(target):
    if target in gene_functions:
        return gene_functions[target]
    return ("Unknown function", "Unknown")

def is_priority(target):
    _, cat = annotate(target)
    return cat in priority_categories

# ─────────────────────────────────────────────
# Analysis per phase
# ─────────────────────────────────────────────
results = {}

for phase, df in [("EXP", exp_named), ("STA", sta_named)]:
    df = df.copy()
    df["function"], df["category"] = zip(*df["target"].apply(annotate))
    df["is_priority"] = df["target"].apply(is_priority)
    df_sorted = df.sort_values("deltaG")

    # Per sRNA: get best (lowest ΔG) interaction, count targets, count priority targets
    srna_summary = []
    for srna, grp in df_sorted.groupby("srna"):
        best_dG = grp["deltaG"].min()
        n_targets = grp["target"].nunique()
        n_priority = grp[grp["is_priority"]]["target"].nunique()
        priority_targets = grp[grp["is_priority"]].sort_values("deltaG")
        srna_summary.append({
            "srna": srna,
            "best_deltaG": best_dG,
            "n_named_targets": n_targets,
            "n_priority_targets": n_priority,
            "priority_targets_df": priority_targets,
            "all_targets_df": grp.sort_values("deltaG"),
        })

    # Sort by: most priority targets, then lowest best_deltaG
    srna_summary.sort(key=lambda x: (-x["n_priority_targets"], x["best_deltaG"]))
    results[phase] = srna_summary

# ─────────────────────────────────────────────
# Print and save summary
# ─────────────────────────────────────────────
output_lines = []

for phase in ["EXP", "STA"]:
    output_lines.append(f"\n{'='*70}")
    output_lines.append(f"  PHASE: {phase}")
    output_lines.append(f"{'='*70}\n")

    for entry in results[phase]:
        srna = entry["srna"]
        output_lines.append(f"sRNA: {srna}")
        output_lines.append(f"  Best ΔG: {entry['best_deltaG']:.2f} kcal/mol")
        output_lines.append(f"  Named DE mRNA targets: {entry['n_named_targets']}")
        output_lines.append(f"  Priority (stress/SMG) targets: {entry['n_priority_targets']}")
        output_lines.append("  All named targets (sorted by ΔG):")
        for _, row in entry["all_targets_df"].iterrows():
            flag = " [PRIORITY]" if row["is_priority"] else ""
            output_lines.append(
                f"    {row['target']:20s}  ΔG={row['deltaG']:7.2f}  | "
                f"{row['function']} | {row['category']}{flag}"
            )
        output_lines.append("")

summary_text = "\n".join(output_lines)
output_file = os.path.join(base_dir, "candidate_summary.txt")

with open(output_file, "w") as f:
    f.write(summary_text)

print(summary_text)
print(f"\n✅ Saved to {output_file}")