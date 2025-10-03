const e=`# Test Article

## 1. Is ionization achieved through direct electron transfer or proton transfer?

**Answer: Both exist, depending on the ionization technique used.**

The essence of mass spectrometry ionization is to impart charge to neutral molecules, turning them into charged ions (positive or negative ions) so that they can move and be detected in magnetic or electric fields. Charge introduction primarily occurs in two ways:

* **Electron transfer (direct gain or loss of electrons):**
  * **Typical example: Electron Impact Ionization Source (EI).** Under vacuum conditions, high-energy electrons (typically 70 eV) emitted from a tungsten filament bombard vaporized sample molecules. After being bombarded, the molecule may directly lose an electron, forming a positively charged molecular ion (M⁺•).
  * \`M + e⁻ → M⁺• + 2e⁻\`
  * This process is high-energy and can easily cause chemical bond cleavage, producing many fragment ions.

* **Proton transfer (gain/loss of protons):**
  * **Typical examples: Electrospray Ionization Source (ESI) and Atmospheric Pressure Chemical Ionization Source (APCI).** These are "soft ionization" techniques, more commonly used in life sciences (e.g., metabolomics, proteomics).
  * In ESI, during the continuous evaporation and Coulombic explosion of charged droplets, **protons (H⁺) are ultimately transferred to the analyte molecule**, forming adduct ions such as \`[M+H]⁺\` (positive ion mode) or are removed from the molecule, forming \`[M-H]⁻\` (negative ion mode).
  * \`M + H⁺ → [M+H]⁺\` (positive ion mode)
  * \`M - H⁺ → [M-H]⁻\` (negative ion mode)
  * Additionally, other cations may be adducted, such as \`[M+Na]⁺\`, \`[M+NH₄]⁺\`, etc.

**Summary:** Traditional EI sources primarily achieve ionization through **electron transfer**; whereas modern soft ionization techniques (e.g., ESI, APCI, MALDI) mainly form ions through **proton transfer** or adduct formation (e.g., with Na⁺).

---

## 2. How to ensure quantitative accuracy?

The instability of ionization degree (i.e., ionization efficiency) is one of the main challenges in mass spectrometry quantitative analysis. The ionization efficiency of the same molecule can vary significantly across different matrices and at different time points. To address this issue, scientists employ the following core strategies to ensure quantitative accuracy:

1. **Internal Standard (IS) Method: This is the most core and effective method.**
    * **Stable Isotope-Labeled Internal Standard (Stable Isotope-Labeled IS):** This is the gold standard. Use a compound with the exact same chemical structure as the analyte but with some atoms replaced by stable isotopes (e.g., ²H, ¹³C, ¹⁵N) as the internal standard (e.g., when detecting amino acids, use ¹³C-labeled versions of the same amino acid). Due to nearly identical chemical properties, the internal standard and analyte behave very similarly during sample preparation and ionization. Any factors causing fluctuations in ionization efficiency affect both the analyte and internal standard equally. By calculating the ratio of **analyte peak area / internal standard peak area**, errors due to ionization efficiency fluctuations can be greatly mitigated.
    * **Structural Analog Internal Standard (Structural Analog IS):** If an isotopic internal standard is unavailable, a compound with a structure very similar to the analyte is chosen as the internal standard, but its correction effect is not as perfect as that of the isotopic internal standard.

2. **Calibration Curve Method:**
    * Prepare a series of standard solutions using known concentrations of standards, adding a **fixed amount of internal standard** to each.
    * After instrument detection, plot a calibration curve (usually a linear regression curve) with the peak area ratio of analyte to internal standard as the y-axis and the standard concentration as the x-axis.
    * The results of unknown samples are calculated using this curve to achieve quantification.

3. **Quality Control (QC) Samples:**
    * Periodically insert QC samples prepared from pooled samples throughout the analytical batch. By monitoring the stability and precision of QC sample results, the reproducibility of the entire experimental process and the reliability of the data can be assessed.

**Summary:** The combination of **stable isotope internal standard + calibration curve method** can maximally correct errors caused by differences in ionization efficiency, sample preparation losses, instrument fluctuations, etc., thereby ensuring quantitative accuracy.

---

## 3. Which ionization sources are used for metabolomics and proteomics, respectively?

These two major omics fields currently heavily rely on **soft ionization techniques** to obtain intact molecular ion information.

* **Metabolomics:**
  * **Primary Ionization Source: Electrospray Ionization Source (ESI).** Due to its excellent "soft" characteristics, it efficiently produces \`[M+H]⁺\` or \`[M-H]⁻\` ions, making it very suitable for analyzing metabolites of varying polarity and small molecular weight (e.g., amino acids, organic acids, sugars, lipids, etc.). It is often coupled with liquid chromatography (LC), i.e., LC-ESI-MS.
  * **Secondary Ionization Sources: Atmospheric Pressure Chemical Ionization Source (APCI)** and **Atmospheric Pressure Photoionization Source (APPI).** APCI has good ionization efficiency for medium to non-polar small molecules (e.g., certain lipids, steroid hormones) and is more tolerant to salts and matrix effects than ESI.

* **Proteomics:**
  * **Primary Ionization Source: Electrospray Ionization Source (ESI).** Especially nano-electrospray (nano-ESI), coupled with high-performance liquid chromatography (particularly nano-flow liquid chromatography, nLC), is the absolute mainstream technology for "bottom-up" proteomics (analysis of enzymatic peptides), known for its extremely high sensitivity.
  * **Important Ionization Source: Matrix-Assisted Laser Desorption/Ionization Source (MALDI).** Commonly used for "peptide mass fingerprinting" identification and imaging mass spectrometry (MSI). It is characterized by simple operation, fast speed, and good salt tolerance, often coupled with time-of-flight mass spectrometry (TOF) (MALDI-TOF/TOF).

---

## 4. What are soft/hard ionization sources?

This classification is based on the amount of energy transferred to the molecule during ionization and whether it easily produces fragments.

* **Hard Ionization Source:**
  * **Characteristics:** The ionization process is high-energy, **easily causing chemical bond cleavage in molecular ions, producing many fragment ions**, while the abundance of the molecular ion itself (parent ion) may be weak or even undetectable.
  * **Advantages:** Fragment ions provide rich structural information, useful for structural elucidation of unknowns and spectral library searching.
  * **Representative: Electron Impact Ionization Source (EI).** The standard ionization source for GC-MS.

* **Soft Ionization Source:**
  * **Characteristics:** The ionization process is gentle, **mainly producing intact molecular ions (or proton adduct ions)**, with few or no fragments.
  * **Advantages:** Very suitable for analyzing complex mixtures, large molecules (e.g., proteins, nucleic acids), and thermally unstable compounds, allowing direct acquisition of accurate molecular weight information.
  * **Disadvantages:** Lack of structural information, usually requiring tandem mass spectrometry (MS/MS) to induce fragmentation for structural information.
  * **Representatives: Electrospray Ionization Source (ESI), Matrix-Assisted Laser Desorption/Ionization Source (MALDI), Atmospheric Pressure Chemical Ionization Source (APCI).**

---

## 5. How to directly obtain the charge number of an ion?

For sources like ESI that can produce multi-charge ions, directly determining the charge number (z) of an ion is key to interpreting its molecular weight (M). There are two main methods:

1. **Isotopic Peak Distribution Method (Isotopic Peak Distribution - the most direct method)**
    * **Principle:** Due to the natural isotopes of elements (e.g., carbon, hydrogen, nitrogen, sulfur, etc.), a molecular ion peak will consist of a cluster of isotopic peaks with a mass difference of 1 Da (e.g., ¹²C and ¹³C). **The mass difference between adjacent isotopic peaks Δm = 1/z.**
    * **Operation:** On a high-mass resolution mass spectrum (e.g., Q-TOF, Orbitrap, FT-ICR), zoom in to observe a peak cluster.
        * Measure the mass difference (Δm) between two adjacent isotopic peaks (e.g., the monoisotopic peak and the M+1 peak).
        * The charge number **z = 1 / Δm**.
    * **Example:** If the measured mass difference between two peaks is 0.5 Da, then the charge number of the ion is z = 1 / 0.5 = **2+**.

2. **Charge Deconvolution Method**
    * **Principle:** The same molecule forms a series of ions with different charge numbers but the same mass (e.g., \`[M+2H]²⁺\`, \`[M+3H]³⁺\`, \`[M+4H]⁴⁺\`...). These ions are distributed on the m/z axis, but through mathematical algorithms, the **neutral molecular weight (M)** they collectively correspond to can be derived. Once M is calculated, combined with the m/z value of any one ion, its charge number z can be determined.
    * **Formula:** For an \`[M+nH]ⁿ⁺\` ion, \`m/z = (M + n*1.007825) / n\`
        * m/z is the measured value
        * M is the neutral molecular weight to be determined
        * n is the charge number (also the number of added protons)
    * **Operation:** Modern mass spectrometer software includes automatic "deconvolution" functions that can automatically identify multiple charge state peaks belonging to the same molecule and directly output its single-charge molecular weight M.

**Summary:** In today's high-resolution mass spectrometry, **directly reading the charge number by analyzing the spacing of isotopic peaks** is the most common and intuitive method.

I hope the above detailed explanations help you fully understand these issues!`;export{e as default};
