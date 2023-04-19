[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![DOI](https://zenodo.org/badge/623305557.svg)](https://zenodo.org/badge/latestdoi/623305557)


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href="https://github.com/FunctionalGeneticsLab/GRADE">
    <!-- <img src="images/pipelineSimplels.png" alt="Logo" width=400>-->
  </a>

  <h3 align="center">GRADE</h3>

  <p align="center">
    General RNAseq Analysis for Differential Expression
    <br />
    <a href="https://github.com/FunctionalGeneticsLab/GRADE"><strong>Explore the docs »</strong></a>
    <br />
    <br />
    <!--<a href="https://github.com/FunctionalGeneticsLab/GRADE">View Demo</a>
    ·--->
    <a href="https://github.com/FunctionalGeneticsLab/GRADE/issues">Report Bug</a>
    ·
    <a href="https://github.com/FunctionalGeneticsLab/GRADE/issues">Request Feature</a>
  </p>
</p>

<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#overview">Overview</a>
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#pipeline-prerequisites">Pipeline Prerequisites</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#license">License</a></li>
  </ol>
</details>

<!-- ABOUT THE PROJECT -->
## Overview

This repository contains the scripts used in our in-house RNAseq analysis pipeline that runs on a PBS cluster.

Each script file foccus on one part of the analysis, from file preparation and quality assessment to differential expression.

0. Script File 0: Prepare pipeline environment.
1. Script File 1: FastQC.
2. Script File 2: Trimmomatic.
3. Script File 3: FastQC.
5. Script File 5: Kallisto Quantification analysis.
6. Script File 6: Kallisto Tables.
8. Script File 8: STAR Alingment.
9. Script File 9: SAM Tools.
10. Script File 10: Novosort.
11. Script File 11: RSEM Quantification Analysis.
12. Script File 12: RSEM Tables.
13. Script File 13: RSEM Transcript to gene.
14. Script File 14: Differential Expression Analysis with EdgeR.

<!-- GETTING STARTED -->
## Pipeline Prerequisites

To get a local copy up and running, make sure you have each script file prerequisites instaled and up to date.

<!--  To get a local copy up and running, make sure you have the following prerequisites instaled and up to date (according to each script file mentioned before):

1. Script File 1 -->

<!-- USAGE EXAMPLES -->
## Usage

Each script can be executed in a PBS cluster by using the following command line:
 
```sh
bash script-name.sh
```

<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to be learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request


<!-- ACKNOWLEDGEMENTS
## Acknowledgements

* []()
* []()
* []() -->


<!-- CONTACT -->
## Contact

Isabela Almeida - mb.isabela42@gmail.com

Project Link: [https://github.com/FunctionalGeneticsLab/GRADE](https://github.com/FunctionalGeneticsLab/GRADE)

<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE` for more information.


<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/FunctionalGeneticsLab/NB-lncRNAs.svg?style=for-the-badge
[contributors-url]: https://github.com/FunctionalGeneticsLab/GRADE/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/FunctionalGeneticsLab/NB-lncRNAs.svg?style=for-the-badge
[forks-url]: https://github.com/FunctionalGeneticsLab/GRADE/network/members
[stars-shield]: https://img.shields.io/github/stars/FunctionalGeneticsLab/NB-lncRNAs.svg?style=for-the-badge
[stars-url]: https://github.com/FunctionalGeneticsLab/GRADE/stargazers
[issues-shield]: https://img.shields.io/github/issues/FunctionalGeneticsLab/NB-lncRNAs.svg?style=for-the-badge
[issues-url]: https://github.com/FunctionalGeneticsLab/GRADE/issues
[license-shield]: https://img.shields.io/github/license/FunctionalGeneticsLab/NB-lncRNAs.svg?style=for-the-badge
[license-url]: https://github.com/FunctionalGeneticsLab/GRADE/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
