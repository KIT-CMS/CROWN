from config import config
import quantities as q

class producer:
    def __init__(
        self,
        call,
        inputs,
        output):
        self.call=call
        self.inputs=inputs
        self.output=output

        #keep track of variable dependencies
        for q in self.inputs:
            q.children.append(self.output)

    def shift(self, name):
        self.output.shift(name)
    
    def writecall(self, shift=""):
        config[shift]["output"]=self.output.get_leaf(shift)
        config[shift]["input"]=",".join([x.get_leaf(shift) for x in self.inputs])
        config[shift]["input_coll"]='{"'+'","'.join([x.get_leaf(shift) for x in self.inputs])+'"}'
        config[shift]["df"]="{df}"
        return self.call.format(**config[shift])

    def writecalls(self):
        calls = [self.writecall()]
        for shift in self.output.shifts:
            calls.append(self.writecall(shift))
        return calls

prod1=producer("phyticsd::FilterID(auto {df}, const std::string {output}, std::string isolationName)", [], q.pt_1)
prod2=producer("FilterID(auto {df}, const std::string {output}, std::string isolationName)", [], q.pt_2)
prod3=producer('FilterID({df}, "{output}", {input_coll}, {ptcut})', [q.pt_1, q.pt_2], q.m_vis)
