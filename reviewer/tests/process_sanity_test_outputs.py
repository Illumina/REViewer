import os
import sys

def main():
    test_output_dir = sys.argv[1]
    frozen_output_dir = sys.argv[2]

    with open('sanity_test_outputs.html', 'w') as html_file:
        html_file.write(
            '''<!DOCTYPE html>
<html>
  <head>
    <title>Sanity test outputs</title>
  </head>
  <style>
  * {
  box-sizing: border-box;
  }
  img{
  max-width: 100%;
  }
  .column {
    float: left;
    position: relative;
    width: 50%;
    height: auto;
    padding: 10px;
  }
  .row {
    border: 1px solid gray;
    height: auto; 
    position: relative;
    overflow: hidden;
  }
  </style>
  <body>''')

        for image_filename in os.listdir(test_output_dir):
            sample_id = image_filename.split(".")[0]
            if not os.path.exists(os.path.join(frozen_output_dir, image_filename)):
                raise ValueError(f"No matching frozen result for test output {image_filename}")
            html_file.write('<div class="row">')

            # New result
            html_file.write('<div class="column">')
            html_file.write(f"<h2>New {sample_id}</h2>")
            html_file.write(f'\n<img src="{os.path.join(test_output_dir, image_filename)}"/>')
            html_file.write('</div>')

            # Frozen result
            html_file.write('<div class="column">')
            html_file.write(f"<h2>Frozen {sample_id}</h2>")
            html_file.write(f'\n<img src="{os.path.join(frozen_output_dir, image_filename)}"/>')
            html_file.write("</div>")

            html_file.write("</div>")
        html_file.write("</body>")
        html_file.write("</html>")
if __name__ == "__main__":
    main()
