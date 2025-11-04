const e=`# Test Article

## 1. Does ionization occur through direct electron transfer or proton transfer?

**Answer: Both exist, depending on the ionization technique used.**

The essence of mass spectrometry ionization is to impart charge to neutral molecules, turning them into charged ions (positive or negative ions) so they can move in magnetic or electric fields and be detected. Charge introduction mainly occurs through two methods:

* **Electron transfer (direct gain/loss of electrons):**
  * **Typical representative: Electron Impact Ion Source (EI).** Under vacuum conditions, high-energy electrons (typically 70 eV) emitted from a tungsten filament bombard vaporized sample molecules. After being bombarded, molecules may directly lose an electron, forming a positively charged molecular ion (M⁺•).
  * \`M + e⁻ → M⁺• + 2e⁻\`
  * This process has high energy and easily leads to chemical bond breakage, producing many fragment ions.

* **Proton transfer (gain/loss of protons):**
  * **Typical representatives: Electrospray Ionization (ESI) and Atmospheric Pressure Chemical Ionization (APCI).** These are "soft ionization" techniques more commonly used in life sciences (such as metabolomics, proteomics).
  * In ESI, during the continuous evaporation and Coulomb explosion of charged droplets, **protons (H⁺) are ultimately transferred to the analyte molecule**, forming adduct ions such as \`[M+H]⁺\` (positive ion mode) or removed from the molecule to form \`[M-H]⁻\` (negative ion mode).
  * \`M + H⁺ → [M+H]⁺\` (positive ion mode)
  * \`M - H⁺ → [M-H]⁻\` (negative ion mode)
  * Additionally, adduction with other cations such as \`[M+Na]⁺\`, \`[M+NH₄]⁺\`, etc., may occur.

**Summary:** Traditional EI sources primarily achieve ionization through **electron transfer**; whereas modern soft ionization techniques (such as ESI, APCI, MALDI) mainly form ions through **proton transfer** or adduct ion formation (e.g., Na⁺).

---

## 2. How to ensure quantitative accuracy?

The instability of ionization efficiency (i.e., ionization efficiency) is the main challenge in mass spectrometry quantitative analysis. The ionization efficiency of the same molecule can vary significantly in different matrices and at different time points. To address this issue, scientists have adopted the following core strategies to ensure quantitative accuracy:

1. **Internal Standard Method (IS): This is the most core and effective method.**
    * **Stable Isotope-Labeled Internal Standard (Stable Isotope-Labeled IS):** This is the gold standard. Use a compound with exactly the same chemical structure as the analyte but with some atoms replaced by stable isotopes (such as ²H, ¹³C, ¹⁵N) as the internal standard (e.g., when detecting amino acids, use ¹³C-labeled same amino acid). Since the chemical properties are almost identical, the behavior of the internal standard and the analyte during sample preparation and ionization is highly consistent. Any factors causing ionization efficiency fluctuations will equally affect the analyte and the internal standard. By calculating the ratio of **analyte peak area / internal standard peak area**, errors caused by ionization efficiency fluctuations can be greatly offset.
    * **Structural Analog Internal Standard (Structural Analog IS):** If no isotope internal standard is available, a compound with a structure very similar to the analyte is chosen as the internal standard, but its correction effect is not as perfect as that of the isotope internal standard.

2. **Calibration Curve Method:**
    * Prepare a series of standard solutions using known concentrations of standards, adding a **fixed amount of internal standard**.
    * After instrument detection, plot a calibration curve (usually a linear regression curve) with the peak area ratio of analyte to internal standard as the y-axis and the standard concentration as the x-axis.
    * The results of unknown samples are calculated through this curve to achieve quantification.

3. **Quality Control (QC) Samples:**
    * Periodically insert QC samples prepared from mixed samples throughout the analysis batch. By monitoring the stability and precision of QC sample results, the repeatability of the entire experimental process and the reliability of the data can be assessed.

**Summary:** Through the combination of **stable isotope internal standard + calibration curve method**, errors caused by differences in ionization efficiency, sample preparation losses, instrument fluctuations, etc., can be corrected to the greatest extent, thereby ensuring quantitative accuracy.

---

## 3. What ionization sources are used for metabolomics and proteomics respectively?

These two major omics fields currently highly rely on **soft ionization techniques** to obtain complete molecular ion information.

* **Metabolomics:**
  * **Main ionization source: Electrospray Ionization (ESI).** Due to its excellent "soft" characteristics, it efficiently produces \`[M+H]⁺\` or \`[M-H]⁻\` ions, making it very suitable for analyzing metabolites with diverse polarities and small molecular weights (such as amino acids, organic acids, sugars, lipids, etc.). It is often coupled with liquid chromatography (LC), i.e., LC-ESI-MS.
  * **Secondary ionization sources: Atmospheric Pressure Chemical Ionization (APCI)** and **Atmospheric Pressure Photoionization (APPI).** APCI has good ionization efficiency for moderately polar and non-polar small molecules (such as certain lipids, steroid hormones) and is more tolerant to salts and matrix effects than ESI.

* **Proteomics:**
  * **Main ionization source: Electrospray Ionization (ESI).** Especially nano-electrospray (nano-ESI), coupled with high-performance liquid chromatography (particularly nano-flow liquid chromatography nLC), is the absolute mainstream technology for "bottom-up" proteomics (enzymatic peptide analysis), characterized by extremely high sensitivity.
  * **Important ionization source: Matrix-Assisted Laser Desorption/Ionization (MALDI).** Commonly used for "peptide mass fingerprint" identification and imaging mass spectrometry (MSI). Its characteristics include simple operation, fast speed, good salt tolerance, and is often coupled with time-of-flight mass spectrometry (TOF) (MALDI-TOF/TOF).

---

## 4. What are soft/hard ionization sources?

This classification is based on the amount of energy delivered to the molecule during ionization and whether it easily produces fragments.

* **Hard Ionization Source:**
  * **Characteristics:** The ionization process has high energy, **easily leading to chemical bond breakage in molecular ions, producing many fragment ions**, while the molecular ion itself (parent ion) may have weak abundance or may not even be detected.
  * **Advantages:** Fragment ions provide rich structural information, which can be used for structural analysis of unknowns and spectral library searching.
  * **Representative: Electron Impact Ion Source (EI).** Standard ionization source for GC-MS.

* **Soft Ionization Source:**
  * **Characteristics:** The ionization process is gentle, **mainly producing intact molecular ions (or proton adduct ions)**, with few or no fragments.
  * **Advantages:** Very suitable for analyzing complex mixtures, large molecules (such as proteins, nucleic acids), and thermally unstable compounds, allowing direct acquisition of accurate molecular weight information.
  * **Disadvantages:** Lack of structural information, usually requiring tandem mass spectrometry (MS/MS) to induce fragmentation for structural information.
  * **Representatives: Electrospray Ionization (ESI), Matrix-Assisted Laser Desorption/Ionization (MALDI), Atmospheric Pressure Chemical Ionization (APCI).**

---

## 5. How to directly obtain the charge number of an ion?

For ionization sources like ESI that can produce multiply charged ions, directly determining the charge number (z) of an ion is key to resolving its molecular weight (M). There are mainly the following two methods:

1. **Isotopic Peak Distribution Method (Isotopic Peak Distribution - Most direct method)**
    * **Principle:** Due to the natural isotopes of elements (such as carbon, hydrogen, nitrogen, sulfur, etc.), a molecular ion peak consists of a cluster of isotopic peaks with a mass difference of 1 Da (e.g., ¹²C and ¹³C). **The mass difference between adjacent isotopic peaks Δm = 1/z.**
    * **Operation:** On a high-mass-resolution mass spectrum (such as Q-TOF, Orbitrap, FT-ICR), zoom in to observe a peak cluster.
        * Measure the mass difference (Δm) between two adjacent isotopic peaks (e.g., the monoisotopic peak and the M+1 peak).
        * Charge number **z = 1 / Δm.**
    * **Example:** If the measured mass difference between two peaks is 0.5 Da, then the charge number of the ion z = 1 / 0.5 = **2+**.

2. **Charge Deconvolution Method (Charge Deconvolution)**
    * **Principle:** The same molecule forms a series of ions with different charge numbers but the same mass (e.g., \`[M+2H]²⁺\`, \`[M+3H]³⁺\`, \`[M+4H]⁴⁺\`...). These ions are distributed on the m/z axis, but through mathematical algorithms, the **neutral molecular weight (M)** they collectively correspond to can be derived. Once M is calculated, combined with the m/z value of any ion, the charge number z can be determined.
    * **Formula:** For \`[M+nH]ⁿ⁺\` ion, \`m/z = (M + n*1.007825) / n\`
        * m/z is the measured value
        * M is the neutral molecular weight to be determined
        * n is the charge number (also the number of protons added)
    * **Operation:** Modern mass spectrometer software comes with automatic "deconvolution" functions that can automatically identify multiple charge state peaks belonging to the same molecule and directly output its single-charge molecular weight M.

**Summary:** In today's high-resolution mass spectrometry, **directly reading the charge number by analyzing the spacing of isotopic peaks** is the most commonly used and intuitive method.

I hope the above detailed answers help you fully understand these questions!`;export{e as default};
