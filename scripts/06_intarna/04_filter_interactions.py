#!/usr/bin/env python3

import os
import pandas as pd

# ─────────────────────────────────────────────
# Paths / settings
# ─────────────────────────────────────────────
base_dir = os.path.join("results", "sRNA_detection", "06_intarna")
out_dir = os.path.join(base_dir, "filtered_interactions")
os.makedirs(out_dir, exist_ok=True)

DG_CUTOFF = -10.0


def find_combined_file(filename):
    """Support combined files either in 06_intarna/ or combined_results/."""
    candidates = [
        os.path.join(base_dir, filename),
        os.path.join(base_dir, "combined_results", filename),
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    raise FileNotFoundError(
        f"Could not find {filename}. Checked:\n  " + "\n  ".join(candidates)
    )


exp_path = find_combined_file("exp_predictions_combined.csv")
sta_path = find_combined_file("sta_predictions_combined.csv")


# ─────────────────────────────────────────────
# Target normalization
# ─────────────────────────────────────────────
def target_to_gene(target):
    """Collapse locus-specific FASTA names, e.g. kfiA__ECOLIN_x -> kfiA."""
    target = str(target).strip()
    if "__ECOLIN_" in target:
        return target.split("__ECOLIN_", 1)[0]
    return target


# ─────────────────────────────────────────────
# Load IntaRNA results
# ─────────────────────────────────────────────
def load_predictions(path):
    df = pd.read_csv(path, sep=";")

    if "E" in df.columns and "deltaG" not in df.columns:
        df.rename(columns={"E": "deltaG"}, inplace=True)

    df = df[df["srna"].astype(str) != "srna"].copy()
    df["deltaG"] = pd.to_numeric(df["deltaG"], errors="coerce")
    df.dropna(subset=["deltaG"], inplace=True)

    df["target"] = df["target"].astype(str).str.strip()
    df = df[~df["target"].str.startswith("N_A")].copy()
    df["target_gene"] = df["target"].apply(target_to_gene)

    return df


exp_df = load_predictions(exp_path)
sta_df = load_predictions(sta_path)

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


# ─────────────────────────────────────────────
# Four broad reporting themes based on Yim et al. (2020)
# ─────────────────────────────────────────────
#
# IMPORTANT:
# Biological relevance is determined by the broader priority_categories above.
# The four Yim themes below are used only to ORGANISE the retained genes.
# Therefore a relevant target such as uvrD is not discarded simply because
# Yim et al. did not explicitly list it under one of their highlighted genes.
#
# Individual gene assignments are based on the existing functional annotation
# category and are broader interpretive grouping, not a claim that Yim et al.
# explicitly assigned every gene to these four groups.

CATEGORY_TO_YIM_GROUP = {
    # Biofilm formation / adhesion / motility / envelope-associated traits
    "Flagellar motility": "Biofilm formation",
    "Chemotaxis / motility": "Biofilm formation",
    "c-di-GMP-mediated motility control": "Biofilm formation",
    "Biofilm / adhesion": "Biofilm formation",
    "Biofilm / curli": "Biofilm formation",
    "Biofilm": "Biofilm formation",
    "Biofilm / c-di-GMP signaling": "Biofilm formation",
    "Virulence / proteolysis": "Biofilm formation",
    "Genotoxin / virulence": "Biofilm formation",
    "Virulence / actin motility": "Biofilm formation",
    "LPS biosynthesis / envelope": "Biofilm formation",
    "Capsule biosynthesis": "Biofilm formation",
    "Envelope stress / peptidoglycan": "Biofilm formation",
    "Protein secretion / envelope": "Biofilm formation",
    "Peptidoglycan remodeling": "Biofilm formation",

    # Stress responses, genome maintenance, protein quality, signalling
    "Stress response": "Stress responses",
    "Heat stress / chaperone": "Stress responses",
    "Acid stress response (GDAR system)": "Stress responses",
    "Acid stress response": "Stress responses",
    "Osmotic stress response": "Stress responses",
    "Oxidative stress response": "Stress responses",
    "Oxidative stress / periplasm": "Stress responses",
    "Protein degradation / stress": "Stress responses",
    "Membrane protein quality control": "Stress responses",
    "Protein quality control / stress": "Stress responses",
    "Heat/oxidative stress": "Stress responses",
    "Membrane stress / stationary phase": "Stress responses",
    "Osmoprotectant synthesis": "Stress responses",
    "Osmoprotectant regulation": "Stress responses",
    "DNA repair / SOS response": "Stress responses",
    "DNA repair": "Stress responses",
    "Ribosome hibernation / stationary phase": "Stress responses",
    "Nucleoid-associated protein / gene regulation": "Stress responses",
    "Two-component signal transduction": "Stress responses",
    "Alkaline stress response": "Stress responses",
    "Cold shock / stress": "Stress responses",

    # Metal acquisition / transport / utilization
    "Iron acquisition": "Metal ion transport/utilization",
    "Zinc homeostasis": "Metal ion transport/utilization",
    "Zinc stress response": "Metal ion transport/utilization",
    "Copper homeostasis": "Metal ion transport/utilization",
    "Copper/silver resistance": "Metal ion transport/utilization",
    "Iron-sulfur cluster assembly": "Metal ion transport/utilization",
    "Iron acquisition (siderophore)": "Metal ion transport/utilization",
    "Iron acquisition / virulence": "Metal ion transport/utilization",
    "Siderophore-mediated virulence": "Metal ion transport/utilization",
    "Siderophore biosynthesis": "Metal ion transport/utilization",

    # Carbon / energy metabolism
    "Anaerobic fermentation": "Carbon metabolism",
    "TCA cycle / anaerobic metabolism": "Carbon metabolism",
    "Anaerobic metabolism": "Carbon metabolism",
    "Respiratory chain": "Carbon metabolism",
    "Aerobic respiration": "Carbon metabolism",
    "TCA cycle": "Carbon metabolism",
    "Hydrogen metabolism": "Carbon metabolism",
    "Formate metabolism": "Carbon metabolism",
    "Anaerobic nitrate respiration": "Carbon metabolism",
    "Fatty acid beta-oxidation": "Carbon metabolism",
}


def annotate(gene):
    return gene_functions.get(gene, ("Unknown function", "Unknown"))


def is_biologically_relevant(gene):
    """Use the original broad biological-relevance criterion."""
    _, category = annotate(gene)
    return category in priority_categories


def yim_group_for_gene(gene):
    """
    Organise a retained relevant gene into one of the four Yim-based themes.
    Returns an empty string only if the detailed category lacks a mapping.
    """
    _, detailed_category = annotate(gene)
    return CATEGORY_TO_YIM_GROUP.get(detailed_category, "")


# Validate the category map so relevant genes are not silently lost.
unmapped_priority_categories = sorted(
    category for category in priority_categories
    if category not in CATEGORY_TO_YIM_GROUP
)

if unmapped_priority_categories:
    raise ValueError(
        "These biologically relevant categories are not mapped to a Yim group:\n  "
        + "\n  ".join(unmapped_priority_categories)
    )


# ─────────────────────────────────────────────
# Filter, classify, and collapse to unique genes
# ─────────────────────────────────────────────
def process_phase(df, phase):
    df = df.copy()

    df["function"], df["detailed_category"] = zip(
        *df["target_gene"].apply(annotate)
    )

    df["passes_energy"] = df["deltaG"] <= DG_CUTOFF
    df["passes_biological_relevance"] = df["target_gene"].apply(
        is_biologically_relevant
    )

    # BOTH criteria define the retained interaction set.
    df["passes_both"] = (
        df["passes_energy"] & df["passes_biological_relevance"]
    )

    # The four Yim themes are applied AFTER retention, for reporting only.
    df["yim_group"] = df["target_gene"].apply(yim_group_for_gene)

    retained_locus = df[df["passes_both"]].copy()

    # Safety check: every retained target must belong to one of the four
    # reporting groups.
    missing_group = retained_locus[retained_locus["yim_group"].eq("")]
    if not missing_group.empty:
        missing = sorted(
            missing_group[["target_gene", "detailed_category"]]
            .drop_duplicates()
            .itertuples(index=False, name=None)
        )
        raise ValueError(
            "Retained genes without a Yim reporting group:\n  "
            + "\n  ".join(f"{g}: {c}" for g, c in missing)
        )

    retained_locus = retained_locus.sort_values(
        ["srna", "yim_group", "deltaG", "target_gene"]
    )

    # Collapse multiple loci with the same gene symbol.
    # Keep the strongest (most negative) ΔG as the representative interaction.
    gene_rows = []
    for (srna, gene), grp in retained_locus.groupby(
        ["srna", "target_gene"], sort=False
    ):
        best = grp.loc[grp["deltaG"].idxmin()]
        loci = sorted(grp["target"].drop_duplicates().astype(str).tolist())

        gene_rows.append({
            "phase": phase,
            "srna": srna,
            "target_gene": gene,
            "best_deltaG": float(best["deltaG"]),
            "yim_group": best["yim_group"],
            "function": best["function"],
            "detailed_category": best["detailed_category"],
            "n_loci": len(loci),
            "target_loci": ", ".join(loci),
        })

    gene_retained = pd.DataFrame(gene_rows)

    if not gene_retained.empty:
        gene_retained = gene_retained.sort_values(
            ["srna", "yim_group", "best_deltaG", "target_gene"]
        )

    # Per-sRNA manuscript summary: counts are UNIQUE TARGET GENES.
    summary_rows = []

    for srna, grp_all in df.groupby("srna"):
        grp_gene = (
            gene_retained[gene_retained["srna"] == srna]
            if not gene_retained.empty
            else pd.DataFrame()
        )

        energy_genes = grp_all.loc[
            grp_all["passes_energy"], "target_gene"
        ].nunique()

        if grp_gene.empty:
            best_prioritized = float("nan")
            targets = ""
        else:
            best_prioritized = grp_gene["best_deltaG"].min()
            targets = ", ".join(
                grp_gene.sort_values("best_deltaG")["target_gene"].tolist()
            )

        summary_rows.append({
            "phase": phase,
            "srna": srna,
            "n_total_target_genes": grp_all["target_gene"].nunique(),
            "n_energy_target_genes_dG_le_minus10": energy_genes,
            "n_prioritized_target_genes": (
                0 if grp_gene.empty
                else grp_gene["target_gene"].nunique()
            ),
            "best_deltaG_overall": grp_all["deltaG"].min(),
            "best_prioritized_deltaG": best_prioritized,
            "prioritized_target_genes": targets,
        })

    summary = pd.DataFrame(summary_rows).sort_values(
        ["n_prioritized_target_genes", "best_prioritized_deltaG"],
        ascending=[False, True],
        na_position="last",
    )

    return retained_locus, gene_retained, summary


exp_locus, exp_gene, exp_summary = process_phase(exp_df, "EXP")
sta_locus, sta_gene, sta_summary = process_phase(sta_df, "STA")


# ─────────────────────────────────────────────
# Save outputs
# ─────────────────────────────────────────────
exp_locus.to_csv(
    os.path.join(out_dir, "exp_filtered_interactions_locus_level.tsv"),
    sep="\t", index=False
)
sta_locus.to_csv(
    os.path.join(out_dir, "sta_filtered_interactions_locus_level.tsv"),
    sep="\t", index=False
)

exp_gene.to_csv(
    os.path.join(out_dir, "exp_filtered_interactions.tsv"),
    sep="\t", index=False
)
sta_gene.to_csv(
    os.path.join(out_dir, "sta_filtered_interactions.tsv"),
    sep="\t", index=False
)

summary = pd.concat([exp_summary, sta_summary], ignore_index=True)
summary.to_csv(
    os.path.join(out_dir, "srna_filtered_interactions_summary.tsv"),
    sep="\t", index=False
)


# ─────────────────────────────────────────────
# Human-readable summary
# ─────────────────────────────────────────────
YIM_GROUP_ORDER = [
    "Biofilm formation",
    "Stress responses",
    "Metal ion transport/utilization",
    "Carbon metabolism",
]

lines = []

for phase, phase_summary, phase_gene in [
    ("EXP", exp_summary, exp_gene),
    ("STA", sta_summary, sta_gene),
]:
    lines.append("=" * 78)
    lines.append(f"PHASE: {phase}")
    lines.append("=" * 78)
    lines.append("")

    for _, row in phase_summary.iterrows():
        srna = row["srna"]

        lines.append(f"sRNA: {srna}")
        lines.append(
            f"  Total phase-matched target genes: "
            f"{int(row['n_total_target_genes'])}"
        )
        lines.append(
            f"  Target genes with ΔG <= -10 kcal/mol: "
            f"{int(row['n_energy_target_genes_dG_le_minus10'])}"
        )
        lines.append(
            f"  Prioritized target genes (energy + biological relevance): "
            f"{int(row['n_prioritized_target_genes'])}"
        )
        lines.append(
            f"  Best ΔG overall: "
            f"{row['best_deltaG_overall']:.2f} kcal/mol"
        )

        if pd.notna(row["best_prioritized_deltaG"]):
            lines.append(
                f"  Best prioritized ΔG: "
                f"{row['best_prioritized_deltaG']:.2f} kcal/mol"
            )

        srna_gene = phase_gene[phase_gene["srna"] == srna]

        for group in YIM_GROUP_ORDER:
            cat = srna_gene[
                srna_gene["yim_group"] == group
            ].sort_values("best_deltaG")

            if cat.empty:
                continue

            lines.append(f"  {group}:")
            for _, hit in cat.iterrows():
                locus_note = (
                    f" [{int(hit['n_loci'])} loci]"
                    if int(hit["n_loci"]) > 1 else ""
                )

                lines.append(
                    f"    {hit['target_gene']:18s} "
                    f"ΔG={hit['best_deltaG']:7.2f} | "
                    f"{hit['function']} | "
                    f"{hit['detailed_category']}"
                    f"{locus_note}"
                )

        lines.append("")

summary_text = "\n".join(lines)

with open(
    os.path.join(out_dir, "candidate_summary_filtered_interactions.txt"),
    "w"
) as handle:
    handle.write(summary_text + "\n")

print(summary_text)
print("=" * 78)
print(f"Saved filtered outputs to: {out_dir}")
print("Counts are unique target genes; duplicate loci are collapsed.")
print("Biological relevance uses the broad priority set.")
print("Yim et al. themes are used for reporting/classification only.")
print("=" * 78)
