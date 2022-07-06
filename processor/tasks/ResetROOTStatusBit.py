import ROOT
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Reset ROOT status bit")
    parser.add_argument("--input", help="input file")
    args = parser.parse_args()
    return args


def reset_status_bit(input_file):
    print(f"Trying to reset status bit for {input_file}")
    rfile = ROOT.TFile(input_file, "UPDATE")
    if "ntuple" not in [x.GetTitle() for x in rfile.GetListOfKeys()]:
        print(f"ntuple tree not found in {input_file}, continueing...")
        rfile.Close()
        return
    t = rfile.Get("ntuple")
    if t.TestBit(ROOT.TTree.EStatusBits.kEntriesReshuffled):
        print("Bit is set, resetting....")
        t.ResetBit(ROOT.TTree.EStatusBits.kEntriesReshuffled)
    rfile.Write()
    rfile.Close()
    print(f"Successfully reset status bit for {input_file}")


# call the function with the input file
if __name__ == "__main__":
    args = parse_args()
    reset_status_bit(args.input)
    print("Done")
    exit(0)
