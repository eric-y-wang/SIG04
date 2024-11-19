# SIG04

### Goal
* Test multiplexed detection of T7oBC5p barcodes
* Test transcriptional activity of library of ligands using Signal-seq
* Test pooled approach (trans signaling)

### Experimental Design
**Cells**: TCRb+CD4+CD44+CD62L-RFP- naive CD4 T cells were sorted from 7W old Foxp3-RFP mice. 2M and 2F mice were harvested. Cells were isolated from skin-draining LNs, mesenteric LNs, and spleen. Cells were activated for 24h with 2.5 ug/mL plate-bound anti-CD3, 3 ug/mL plate-bound anti-CD28, and 30 U/mL hIL2. After 22h, cells were lifted from plate and used in assay.

**Treatment**: Cells were treated in 96-well plates at 100,000 cells per well in 100 ul cRPMI. Cells were spinfected with virus or recombinant proteins 2000g x 20m and viral media was removed. Individually incubated cells were resuspended in 200 uL cRPMI + 30 U/mL hIL2 and transferred to new 96-well plates with anti-CD3/28. Pooled cells were processed with pooling protocol and then transferred to new plates. Cells were incubated in TC incubator for 6.5h. Afterwards, cells were pooled for scRNA-seq.

**Sequencing**: GEM-X 5' scRNA-seq was performed. Barcodes were amplified by custom protocol.
