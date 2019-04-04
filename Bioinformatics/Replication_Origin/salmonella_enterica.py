import urllib.request
import approximate_patterns
import skew_diagram

if __name__ == "__main__":
    genome = urllib.request.urlopen("http://bioinformaticsalgorithms.com/data/Salmonella_enterica.txt").read()
    # decode from byte to string
    genome = genome.decode("utf-8")
    # convert string to list, also removing new lines
    genome = genome.split("\r\n")
    # convert list back to string, without the first and last lines
    genome = "".join(genome[1:len(genome)-1])
    dnaa_boxes = []
    for potential_ori in skew_diagram.minimum_skew(genome):
        ori_window = genome[(potential_ori - 500):(potential_ori + 500)]
        k = 9
        d = 1
        d_box = approximate_patterns.most_frequent_approx_pattern(ori_window, k, d)
        dnaa_boxes.append(d_box)
    print(dnaa_boxes)
