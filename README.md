# LCMS for Galaxy #

Galaxy tool for LC-MS data deisotoping and annotation for metabolomics
analysis.


The original codes by Zoe Hall is located in
[lcms_processing](https://github.com/hallz/lcms_processing_).

## Installation ##

- Install [Galaxy](https://github.com/galaxyproject/galaxy) under Linux.

- Use `git` to clone this tool

  ```bash
  git clone https://github.com/wanchanglin/lcms.git
  ```

- Add this tool's location into Galaxy' tool config file:
  `~/Galaxy/config/tool_conf.xml`. For example, one simplified
  `tool_conf.xml` looks like:

  ```xml
  <?xml version='1.0' encoding='utf-8'?>
  <toolbox monitor="true">

    <section id="getext" name="Get Data">
      <tool file="data_source/upload.xml" />
    </section>

    <section id="MyTools" name="My Tools">
      <tool file="/path/to/lcms/lcms.xml" />
    </section>

  </toolbox>
  ```

- Test data are in `test-data`.

## Authors, contributors & contacts ##

- Wanchang Lin (wl361@cam.ac.uk), University of Cambridge
- Zoe Hall (zlh22@cam.ac.uk), University of Cambridge
- Julian L Griffin (jlg40@cam.ac.uk), University of Cambridge

