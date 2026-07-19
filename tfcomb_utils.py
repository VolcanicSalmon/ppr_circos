class Trigger_Target(list):
    bwpath=None
    def as_tab(self):
        table=[]
        for p in self:
            attrs=getattr(p)
            site1={f""}
