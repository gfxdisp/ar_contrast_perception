# AR Contrast Perception
### [Paper](https://dongyeon93.github.io/assets/pdf/AR_contrast.pdf) | [Project Page](https://www.cl.cam.ac.uk/research/rainbow/projects/ar-contrast/) | [Video](https://www.youtube.com/watch?v=tPUkpA6VI8c)

[Dongyeon Kim](https://dongyeon93.github.io/),
[Maliha Ashraf](https://malihaashraf.github.io/),
[Alexandre Chapiro](https://achapiro.github.io/),
and [Rafał K. Mantiuk](https://www.cl.cam.ac.uk/~rkm38/)

Test code & data for the SIGGRAPH Asia 2025 conference proceeding titled "Supra-threshold Contrast Perception in Augmented Reality"

## 1. To Start

To check out
```
git clone https://github.com/gfxdisp/ar_contrast_perception.git
cd ar_contrast_perception
```

## 2-1. (Optional) Add the CastleCSF Module
Add the castleCSF module as a Git submodule:

```bash
git submodule add https://github.com/gfxdisp/castleCSF.git
```

## 2-2. Submodule update
```
git submodule update --init --recursive
```

## 3. Test & Plot

This is a simple MATLAB example demonstrating how to generate some of the figures in the paper and supplementary.

### 🚀 How to Run

1. Open MATLAB.
2. Navigate to the folder using the command window or `cd`:

    ```matlab
    cd path/to/this/folder
    ```

3. Run the script (Exp1-result, Fig. 5):

    ```matlab
    plot_exp1_main
    ```

4. Run the script (Exp2-model, Fig. 6):

    ```matlab
    plot_ar_contrast_models
    ```

5. Run the script (Exp2-result, Fig. 7):

    ```matlab
    plot_exp2_main
    ```

6. Run the script (Exp3-result, Fig. 8 & Fig.S8):

    ```matlab
    plot_exp3_main
    plot_exp3_per_image
    ```

7. Run the script (Auto-brightness curve, Fig. 9):

    ```matlab
    plot_autobrightness
    ```
---

## 4. 📁 Data

- `data/exp1/exp1_average_data.csv`: Processed data with all observers for Experiment 1
- `data/exp1/exp1_individual_data.csv`: Processed data with individual observers for Experiment 1
- `data/exp2/exp2_average_data.csv`: Processed data with all observers for Experiment 2
- `data/exp2/exp2_individual_data.csv`: Processed data with individual observers for Experiment 2
- `data/exp3/exp3_raw_data.csv`: Raw data with all observers for Experiment 3
- `image/{ImageName}.png`: Images used for Experiment 3

### FYI
#### {observer}
The ids of all observers were anonymized per experiment. Thus, o1 in exp1 is not equal to o1 in exp2.

### Experiment 1
#### {bkg_type}
- `VR-flat-normal-bkg-same-bkg-contrast-1-opposite-static`: SP
- `flat-normal-bkg-same-bkg-contrast-1-opposite-static`: DF
- `optics-normal-bkg-same-bkg-contrast-1-opposite-static`: DNS
- `noise-normal-bkg-same-bkg-contrast-1-opposite-static`: DNP
- `optics-normal-bkg-same-bkg-contrast-1-opposite-slow-dynamic`: DND

### Experiment 2
#### {col_direction}
- `1`: Ach, `2`: RG, `3`: YV

#### {bkg_type}
- `max-flat-normal-bkg-same-bkg-contrast-1-opposite-static`: Ach: c=0.8 / RG: c=0.1 / YV: c=0.8
- `margin-flat-normal-bkg-same-bkg-contrast-1-opposite-static`: Ach: c=0.3 / RG: c=0.07 / YV: c=0.3

## Citation
To be updated

## Contact
If you have any questions, please contact

- Dongyeon Kim (dk721@cam.ac.uk)
- Maliha Ashraf (ma905@cam.ac.uk)
- Rafal Mantiuk (rafal.mantiuk@cl.cam.ac.uk)
