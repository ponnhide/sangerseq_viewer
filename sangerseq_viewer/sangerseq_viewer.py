import os 
import sys 
import copy
import argparse
import patchworklib as pw
from QUEEN.queen import * 
from Bio import SeqIO
from Bio import pairwise2
import matplotlib
matplotlib.rcParams["xtick.major.width"] = 0.4
matplotlib.rcParams["xtick.minor.width"] = 0.4
matplotlib.rcParams["ytick.major.width"] = 0.4
matplotlib.rcParams["ytick.minor.width"] = 0.4

margin     = pw.param["margin"]
_atgc_dict = {0:"A", 1:"T", 2:"G", 3:"C"}
color_dict = {"G":"#f2f059", "C":"#74b2d7", "A":"#79E5B7", "T":"#ff776c", "N":"#FFFFFF", "-":"#FFFFFF", "/":"#FFFFFF"}

def colorbar(ax, color_dict, query, subject, char=False, fontsize=10, label=False, zero_position=0):
    if len(query) > 200:
        bars = ax.bar(list(range(len(query))), [1.0] * (len(query)), width=1.0, edgecolor="#BBBBBB", linewidth=0.1, align="edge",bottom=0.05)
    else:
        bars = ax.bar(list(range(len(query))), [1.0] * (len(query)), width=1.0, edgecolor="#BBBBBB", linewidth=0.1, align="edge",bottom=0.05)
    ax.set_xlim(0,len(query))
    ax.set_ylim(0,1.00)
    p = 0
    for bar, q, s in zip(bars, query, subject):
        color = color_dict[q]
        if q != s:
            #bar.set_edgecolor("#E83929")
            if char == True:
                if q == "/":
                    bar.set_facecolor("#BBBBBB")
                    ax.text(p+0.5,0.45,"-",va="center",ha="center",fontsize=fontsize,zorder=100,color="k")     
                else:
                    bar.set_facecolor("#993366")
                    bar.set_alpha(0.2)
                    if q == "-":
                        ax.text(p+0.5,0.45,"-",va="center",ha="center",fontsize=fontsize,zorder=100,color="k")#fontweight="bold")
                    else:
                        ax.text(p+0.5,0.45,q,va="center",ha="center",fontsize=fontsize,zorder=100,color="k")#fontweight="bold")
        else:
            bar.set_facecolor("#FFFFFF")
            if char == True:
                ax.text(p+0.5,0.45,q,va="center",ha="center",fontsize=fontsize,zorder=100)     
        p += 1 
    
    #ax.set_xticks([])
    if label == True:
        positions  = [pos - 0.5 for pos in range(0, len(query)+1, 10)]
        ticklabels = [str(int(zero_position+pos+0.5)) for pos in positions] 
        if positions[0] < 1 and zero_position == 0:
            positions[0]  = 0.5
            ticklabels[0] = 1
        ax.set_xticks(positions)
        ax.set_xticklabels(ticklabels)
    else:
        ax.set_xticks([])
    
    ax.set_zorder(2.0)
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_linewidth(0.4)  
    #ax.spines["top"].set_visible(False)
    ax.patch.set_alpha(1.0)
    ax.set_xlim(0, len(query))
    return bars

def abi_to_dict(filename):
    record   = SeqIO.read(filename,'abi')
    abi_data = {"conf":[[],[]],
                "channel": {
                           "A":[[],[]],
                           "T":[[],[]],
                           "G":[[],[]],
                           "C":[[],[]],
                           },
                "_channel":{
                            "A":[[],[]],
                            "T":[[],[]],
                            "G":[[],[]],
                            "C":[[],[]], 
                           }
                }
    
    pre_pos = 0 
    pos_set = [] 
    for i, (pos, conf) in enumerate(zip(record.annotations['abif_raw']["PLOC1"], record.annotations['abif_raw']["PCON1"])):
        abi_data["conf"][0].append(i)
        abi_data["channel"]["G"][0].append(i)
        abi_data["channel"]["A"][0].append(i)
        abi_data["channel"]["T"][0].append(i)
        abi_data["channel"]["C"][0].append(i) 
        
        abi_data["conf"][1].append(conf)
        abi_data["channel"]["G"][1].append(record.annotations['abif_raw']["DATA9"][pos])
        abi_data["channel"]["A"][1].append(record.annotations['abif_raw']["DATA10"][pos])
        abi_data["channel"]["T"][1].append(record.annotations['abif_raw']["DATA11"][pos])
        abi_data["channel"]["C"][1].append(record.annotations['abif_raw']["DATA12"][pos])    
       
        step = 0.1
        for j in range(10):
            abi_data["_channel"]["G"][0].append(i+step*j)   
            abi_data["_channel"]["A"][0].append(i+step*j)
            abi_data["_channel"]["T"][0].append(i+step*j) 
            abi_data["_channel"]["C"][0].append(i+step*j)   
            try:
                abi_data["_channel"]["G"][1].append(record.annotations['abif_raw']["DATA9"][pos-5+j])   
                abi_data["_channel"]["A"][1].append(record.annotations['abif_raw']["DATA10"][pos-5+j])
                abi_data["_channel"]["T"][1].append(record.annotations['abif_raw']["DATA11"][pos-5+j]) 
                abi_data["_channel"]["C"][1].append(record.annotations['abif_raw']["DATA12"][pos-5+j])  
            except:
                abi_data["_channel"]["G"][1].append(0)   
                abi_data["_channel"]["A"][1].append(0)
                abi_data["_channel"]["T"][1].append(0) 
                abi_data["_channel"]["C"][1].append(0)  

        pos_set.append((pre_pos, pos))
        pre_pos = pos
    pos_set.append((pre_pos, record.annotations['abif_raw']["PLOC1"][-1]))
    
    return abi_data 

def generate_consensusseq(abidata):
    consensus_seq = ""  
    for values in zip(abidata["channel"]["A"][1], abidata["channel"]["T"][1], abidata["channel"]["G"][1], abidata["channel"]["C"][1]):
        consensus_seq += _atgc_dict[values.index(max(values))]
     
    return (consensus_seq, consensus_seq.translate(str.maketrans("ATGC","TACG"))[::-1]) 

def gl_alignment(template, query, single=False): 
    subject   = template.seq
    alignment = pairwise2.align.globalms(subject+subject, query, 2, 0, -10, -1, penalize_end_gaps=False, one_alignment_only=True)[0] 
    seqB  = alignment.seqB
    seqA  = alignment.seqA
    score = alignment.score
    for s, q in enumerate(seqB):
        if q == "-" or q == "/":
            pass 
        else:
            break 
    start = s
    
    for e, q in enumerate(seqB[::-1]):
        if q == "-" or q == "/":
            pass 
        else:
            break 
    
    end = len(2*subject) - e
    if start > len(subject):
        new_start = start - len(subject)
    else:
        new_start = start
    
    if end > len(subject): 
        new_end = end - len(subject) 
    else:
        new_end = end  

    if single == True:
        if new_start > new_end:
            nongap = joindna(cropdna(template, new_start, len(subject), quinable=0), cropdna(template, 0, new_end, quinable=0), quinable=0)
        else:
            nongap = cropdna(template, new_start, new_end) 
        
        op = 0
        flag = 0 
        gap_positions = []
        ops           = []
        for p, s in enumerate(seqA[start:-1*e]):
            if s == "-" or s == "/":
                if flag == 0:
                    gapstart    = p
                    gapstart_op = op
                    flag = 1
                else:
                    pass 

            else:
                if flag == 1:
                    gapend = p
                    gap_positions.append((gapstart, gapend))
                    ops.append(gapstart_op) 
                    flag = 0
                else:
                    pass 
                op += 1

        if len(gap_positions) == 0:
            wgap = nongap 
        else: 
            pre_op = 0
            for p, (op, gap_pos) in enumerate(zip(ops,gap_positions)):
                if p == 0:
                    wgap = joindna(cropdna(nongap, 0, op), QUEEN(seq="-"*(gap_pos[1]-gap_pos[0])))    
                else:
                    wgap = joindna(wgap, cropdna(nongap, pre_op, op), QUEEN(seq="-"*(gap_pos[1]-gap_pos[0])))    
                pre_op = op
            
            if pre_op - len(nongap.seq) == 0:
                pass
            else:
                wgap = joindna(wgap, cropdna(nongap, pre_op, len(nongap.seq)))
        
        return wgap, seqB[start:-1*e], nongap, new_start, new_end, alignment.score
    
    else:
        return seqB[start:-1*e], new_start, new_end, alignment.score    

def recursive_alignment(template, query_list):
    starts         = [] 
    ends           = []  
    new_query_list = []
    strands        = [] 
    for query in query_list: 
        result_f = gl_alignment(template, query[0], single=False) 
        result_r = gl_alignment(template, query[1], single=False)
        if result_f[-1] >= result_r[-1]:
            aligned_q, start, end, score = result_f
            strands.append(1) 
        else:
            aligned_q, start, end, score = result_r
            strands.append(-1) 
        
        #aligned_q, start, end, score = alignment(template, query, single=False)
        new_query_list.append(aligned_q) 
        starts.append(start) 
        ends.append(end) 
         
    start = min(starts)
    end   = max(ends) 
    if start > end:
        nongap = joindna(cropdna(template, start, len(template), quinable=0), cropdna(template, 0, end, quinable=0), quinable=0)
    else:
        nongap = cropdna(template, start, end)
    
    new_temp = nongap.seq
    while 1:
        templates          = []
        new_new_query_list = [] 
        for query in new_query_list:
            alignment = pairwise2.align.globalms(new_temp, query, 2, 0, -10, -1, penalize_end_gaps=False, one_alignment_only=True)[0]
            templates.append(alignment.seqA) 
            new_new_query_list.append(alignment.seqB)  
        
        new_query_list = new_new_query_list 
        if len(set(templates)) == 1:
            new_temp = templates[0]
            break
        else:
            templates.sort(key=lambda x: len(x))
            new_temp = templates[-1]
    
    op = 0
    flag = 0 
    gap_positions = []
    ops           = []
    for p, s in enumerate(new_temp):
        if s == "-" or s == "/":
            if flag == 0:
                gapstart    = p
                gapstart_op = op
                flag = 1
            else:
                pass 

        else:
            if flag == 1:
                gapend = p
                gap_positions.append((gapstart, gapend))
                ops.append(gapstart_op) 
                flag = 0
            else:
                pass 
            op += 1

    if len(gap_positions) == 0:
        wgap = nongap 
    else: 
        pre_op = 0
        for p, (op, gap_pos) in enumerate(zip(ops,gap_positions)):
            if p == 0:
                wgap = joindna(cropdna(nongap, 0, op), QUEEN(seq="-"*(gap_pos[1]-gap_pos[0])))    
            else:
                wgap = joindna(wgap, cropdna(nongap, pre_op, op), QUEEN(seq="-"*(gap_pos[1]-gap_pos[0])))    
            pre_op = op
        
        if pre_op - len(nongap.seq) == 0:
            pass
        else:
            wgap = joindna(wgap, cropdna(nongap, pre_op, len(nongap.seq)))
    
        
    return wgap, list(zip(new_query_list, strands)), nongap, start, end

def reform_abidata(abi_data, seqB, strand):
    if strand == -1:
        new_channel = copy.deepcopy(abi_data["_channel"]) 
        new_channel["G"][0] = abi_data["_channel"]["C"][0] 
        new_channel["G"][1] = abi_data["_channel"]["C"][1][::-1]
        new_channel["A"][0] = abi_data["_channel"]["T"][0] 
        new_channel["A"][1] = abi_data["_channel"]["T"][1][::-1]
        new_channel["T"][0] = abi_data["_channel"]["A"][0] 
        new_channel["T"][1] = abi_data["_channel"]["A"][1][::-1]
        new_channel["C"][0] = abi_data["_channel"]["G"][0] 
        new_channel["C"][1] = abi_data["_channel"]["G"][1][::-1]
        abi_data["_channel"] = new_channel
        
        new_channel = copy.deepcopy(abi_data["channel"]) 
        new_channel["G"][0] = abi_data["channel"]["C"][0] 
        new_channel["G"][1] = abi_data["channel"]["C"][1][::-1]
        new_channel["A"][0] = abi_data["channel"]["T"][0] 
        new_channel["A"][1] = abi_data["channel"]["T"][1][::-1]
        new_channel["T"][0] = abi_data["channel"]["A"][0] 
        new_channel["T"][1] = abi_data["channel"]["A"][1][::-1]
        new_channel["C"][0] = abi_data["channel"]["G"][0] 
        new_channel["C"][1] = abi_data["channel"]["G"][1][::-1]
        abi_data["channel"] = new_channel
        
        abi_data["conf"][0] = abi_data["conf"][0][::-1]
        abi_data["conf"][1] = abi_data["conf"][1][::-1]

    op   = 0 
    flag = 0 
    gap_positions = []
    ops           = []
    for p, s in enumerate(seqB):
        if s == "-" or s == "/":
            if flag == 0:
                gapstart    = p
                gapstart_op = op
                flag = 1
            else:
                pass 

        else:
            if flag == 1:
                gapend = p
                gap_positions.append((gapstart, gapend)) 
                ops.append(gapstart_op) 
                flag = 0
            else:
                pass 
            op += 1
    #print(seqB)  
    #rrint(gap_positions) 
    abi_data["_channel_wgap"] = {
                                 "A":[[],[]], 
                                 "T":[[],[]],
                                 "G":[[],[]],
                                 "C":[[],[]]
                                }
    opi    = 0 
    gapsum = 0 
    flag   = 0
    if len(ops) > 0:
        op     = ops[opi] * 10
        for i, v in enumerate(abi_data["_channel"]["G"][0]):
            #if i % 10 == 0 or i < 10:
            #    print(i, gapsum) 
            abi_data["_channel_wgap"]["G"][0].append(abi_data["_channel"]["G"][0][i] + gapsum) 
            abi_data["_channel_wgap"]["A"][0].append(abi_data["_channel"]["A"][0][i] + gapsum) 
            abi_data["_channel_wgap"]["T"][0].append(abi_data["_channel"]["T"][0][i] + gapsum)
            abi_data["_channel_wgap"]["C"][0].append(abi_data["_channel"]["C"][0][i] + gapsum) 
            
            abi_data["_channel_wgap"]["G"][1].append(abi_data["_channel"]["G"][1][i]) 
            abi_data["_channel_wgap"]["A"][1].append(abi_data["_channel"]["A"][1][i]) 
            abi_data["_channel_wgap"]["T"][1].append(abi_data["_channel"]["T"][1][i])
            abi_data["_channel_wgap"]["C"][1].append(abi_data["_channel"]["C"][1][i]) 

            if i >= op and flag == 0:
                gaplength = gap_positions[opi][1] - gap_positions[opi][0]
                for j in range(gaplength*10):
                    abi_data["_channel_wgap"]["G"][0].append(abi_data["_channel"]["G"][0][i] + gapsum + j*0.1) 
                    abi_data["_channel_wgap"]["A"][0].append(abi_data["_channel"]["A"][0][i] + gapsum + j*0.1) 
                    abi_data["_channel_wgap"]["T"][0].append(abi_data["_channel"]["T"][0][i] + gapsum + j*0.1)
                    abi_data["_channel_wgap"]["C"][0].append(abi_data["_channel"]["C"][0][i] + gapsum + j*0.1) 
                    
                    abi_data["_channel_wgap"]["G"][1].append(0) 
                    abi_data["_channel_wgap"]["A"][1].append(0) 
                    abi_data["_channel_wgap"]["T"][1].append(0)
                    abi_data["_channel_wgap"]["C"][1].append(0) 
                   
                gapsum += gaplength 
                opi    += 1
                if opi >= len(ops):
                    flag = 1
                    pass
                else:
                    op = ops[opi] * 10
    else:
        abi_data["_channel_wgap"] = copy.deepcopy(abi_data["_channel"]) 
    return abi_data         

def visualize_abi(pos, subject, abi_data, query, abiname=None, zero_position=0, label=False, single=False):
    values = [] 
    ax  = pw.Brick(figsize=(pos.x1 - pos.x0, 0.6)) 
    axpos = ax.get_position() 
    
    ax2 = pw.Brick(figsize=(pos.x1 - pos.x0, 0.6))
    ax2.set_position([axpos.x0, axpos.y0, abs(axpos.x0-axpos.x1), abs(axpos.y0-axpos.y1)]) 
    
    ax.set_zorder(1) 
    ax2.set_zorder(0) 
    
    ax.plot(abi_data["_channel_wgap"]["G"][0], abi_data["_channel_wgap"]["G"][1], color="#e3e14f", lw=1, zorder=0) 
    ax.plot(abi_data["_channel_wgap"]["A"][0], abi_data["_channel_wgap"]["A"][1], color="#79E5B7", lw=1, zorder=0)
    ax.plot(abi_data["_channel_wgap"]["T"][0], abi_data["_channel_wgap"]["T"][1], color="#ff776c", lw=1, zorder=0) 
    ax.plot(abi_data["_channel_wgap"]["C"][0], abi_data["_channel_wgap"]["C"][1], color="#74b2d7", lw=1, zorder=0)  
    values.extend(abi_data["_channel_wgap"]["G"][1]) 
    values.extend(abi_data["_channel_wgap"]["A"][1]) 
    values.extend(abi_data["_channel_wgap"]["T"][1])
    values.extend(abi_data["_channel_wgap"]["C"][1])
    
    try:
        top    = abs(max(values)-min(values))*1.02
        bottom = -1.0*abs(max(values)-min(values))*0.02
    except:
        top    = 1000
        bottom = -20

    i = 0 
    for p, (q,s) in enumerate(zip(query, subject.seq)):
        if q != "/" and q != s:
            ax.bar([p+0.5], [top-bottom], bottom=bottom, width=1, lw=0.0, facecolor="#993366", edgecolor="#DDDDDD", zorder=1, alpha=0.2)

        if q != "-" and q != "/":
            ax2.bar([p+0.5], [abi_data["conf"][1][i]], bottom=0, width=1, lw=0.1, facecolor="#EFEFFF", edgecolor="#DDDDDD", zorder=1)
            i += 1

    ax.spines["top"].set_visible(False)
    ax.spines["bottom"].set_visible(False)  
    ax.spines["right"].set_visible(False) 
    ax.spines["left"].set_visible(False) 
    ax.patch.set_alpha(0.0) 
    if abiname is None:
        pass
    else:
        ax.yaxis.set_label_position("right")
        ax.set_ylabel(abiname, rotation=0, labelpad=10, ha="left")
    
    ax2.spines["top"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)  
    ax2.spines["right"].set_visible(False) 
    ax2.spines["top"].set_linewidth(0.4)  
    ax2.spines["bottom"].set_linewidth(0.4)  
    ax2.spines["right"].set_linewidth(0.4)  
    ax2.spines["left"].set_linewidth(0.4)  
    #ax2.patch.set_alpha(0.0) 

    ax.set_xlim(0, len(subject.seq)) 
    ax.set_xticks([]) 
    ax.set_ylim(bottom, top) 
    ax.set_yticks([])
    
    #ax2.plot([0,len(subject.seq)], [20,20], lw=0.2, ls="--", color="#BBBBBB")
    ax2.set_xlim(0, len(subject.seq))
    ax2.set_xticks([])
    ax2.set_ylim(0,70)
    ax2.set_yticks([0,30,60])
    if single == True:
        ax2.set_ylabel("Quality")
    axquery = pw.Brick(figsize=(pos.x1 - pos.x0, pos.y1 - pos.y0)) 
    colorbar(axquery, color_dict, query, subject.seq, char=True, fontsize=10, label=label, zero_position=zero_position)
    return ax, ax2, axquery 

def visualize(subject, abi_data, query, abiname=None, start=0, end=None, display_quality=True):
    if start == 0 and end is None:
        axmap  = visualizemap(subject, seq=True, linebreak=len(subject.seq)+1, tick_interval=10, title="", height_scale=0.8, fontsize_nucl=10, fontsize=10) 
    else: 
        axmap   = visualizemap(subject, seq=True, start=start, end=end, linebreak=end-start, tick_interval=10, title="", height_scale=0.8, fontsize_nucl=10, fontsize=10) 
        subject = cropdna(subject, start=start, end=end, quinable=0)
        if type(query) in (list, tuple):
            new_query_list = []
            abidata_list   = [] 
            for aquery, aabi_data in zip(query, abi_data):
                conf0 = []
                conf1 = [] 
                i = 0 
                for q in aquery[:start]:
                    if q != "-" and q != "/":
                        i += 1
                    
                j = 0 
                for q in aquery[start:end]:
                    if q != "-" and q != "/":
                        conf0.append(aabi_data["conf"][0][i+j])
                        conf1.append(aabi_data["conf"][1][i+j])
                        j += 1


                new_abidata = {"conf":[conf0, conf1], "_channel_wgap":{}}
                for nucl in "ATGC":
                    positions = [p-start for p in aabi_data["_channel_wgap"][nucl][0] if start <= p <= end]
                    values    = [v for p, v in zip(aabi_data["_channel_wgap"][nucl][0], aabi_data["_channel_wgap"][nucl][1]) if start <= p <= end]
                    new_abidata["_channel_wgap"][nucl] = [positions, values]
                new_query_list.append(aquery[start:end]) 
                abidata_list.append(new_abidata)
            query    = new_query_list
            abi_data = abidata_list
        
        else:
            conf0 = []
            conf1 = [] 
            i = 0 
            for q in query[:start]:
                if q != "-" and q != "/":
                    #conf0.append(abi_data["conf"][0][i+j])
                    #conf1.append(abi_data["conf"][1][i+j])
                    i += 1
                
            j = 0 
            for q in query[start:end]:
                if q != "-" and q != "/":
                    conf0.append(abi_data["conf"][0][i+j])
                    conf1.append(abi_data["conf"][1][i+j])
                    j += 1

            new_abidata = {"conf":[conf0, conf1], "_channel_wgap":{}}
            for nucl in "ATGC":
                positions = [p-start for p in abi_data["_channel_wgap"][nucl][0] if start <= p <= end]
                values    = [v for p, v in zip(abi_data["_channel_wgap"][nucl][0], abi_data["_channel_wgap"][nucl][1]) if start <= p <= end]
                new_abidata["_channel_wgap"][nucl] = [positions, values]
            abi_data = new_abidata
            query    = query[start:end] 

    keys   = list(axmap.bricks_dict.keys())
    axmap.bricks_dict[keys[1]].set_xticks([])
 
    pos     = axmap.bricks_dict[keys[1]].get_position() 
    patches = axmap.bricks_dict[keys[1]].containers[0].patches
    for s, patch in zip(subject.seq, patches):
        patch.set_alpha(0.6) 
        patch.set_facecolor(color_dict[s]) 
        #patch.set_facecolor("k")
 
    if type(query) in (list, tuple):
        pw.param["margin"] = None
        mappos = axmap.bricks_dict[keys[0]].get_position()
        axmap.bricks_dict[keys[0]].change_plotsize([mappos.x1-mappos.x0, (mappos.y1-mappos.y0)*0.95])
        axmap = axmap.bricks_dict[keys[0]]/axmap.bricks_dict[keys[1]]

        subax_list = []  
        for i, (aabi_data, aquery) in enumerate(zip(abi_data, query)):
            if i == len(query) - 1:
                ax1, ax2, axquery  = visualize_abi(pos, subject, aabi_data, aquery, abiname=abiname[i], zero_position=start, label=True) 
            else:
                ax1, ax2, axquery  = visualize_abi(pos, subject, aabi_data, aquery, abiname=abiname[i], zero_position=start, label=False)
            
            if display_quality == True:
                ax = pw.Bricks({ax1.get_label():ax1, ax2.get_label():ax2})
            else:
                ax1.patch.set_alpha(1.0) 
                ax = ax1
            
            subax = (ax/axquery) 
            subax_list.append(subax) 
        pw.param["margin"] = 0.05
        sub_axes = pw.stack(subax_list, operator="/")
        sub_axes.set_supylabel("Quality", labelpad=4) 
        sub_axes.set_supspine()
        sub_axes._case.spines["left"].set_linewidth(0.4) 
        pw.param["margin"] = None
        ax_all = axmap/sub_axes
        pw.param["margin"] = 0.05
    else:
        pw.param["margin"] = None
        mappos = axmap.bricks_dict[keys[0]].get_position()
        axmap.bricks_dict[keys[0]].change_plotsize([mappos.x1-mappos.x0, (mappos.y1-mappos.y0)*0.95])
        axmap = axmap.bricks_dict[keys[0]]/axmap.bricks_dict[keys[1]]
        ax1, ax2, axquery = visualize_abi(pos, subject, abi_data, query, abiname=abiname, zero_position=start, label=True, single=True) 
        
        if display_quality == True:
            ax = pw.Bricks({ax1.get_label():ax1, ax2.get_label():ax2})
        else:
            ax1.patch.set_alpha(1.0) 
            ax = ax1
        ax_all = axmap/(ax/axquery) 
    return ax_all

def view_sanger(gbkpath, abipath, start=None, end=None, linebreak=None, output=None, display_quality=True, dpi=None):
    template = QUEEN(record=gbkpath)
    project  = template.project    
    if start is None and end is None:
        pass 
    else:
        if start is None:
            start = 0    
        if end is None:
            end = len(template.seq) 
    
    if os.path.isdir(abipath):
        abidata_list = [] 
        query_list   = [] 
        abifile_path_list = [] 
        for abifile_path in os.listdir(abipath):
            if ".ab1" == abifile_path[-4:]:
                abifile_path = abipath + "/" + abifile_path 
                abidata      = abi_to_dict(abifile_path)
                query        = generate_consensusseq(abidata)
                abidata_list.append(abidata)             
                query_list.append(query)
                abifile_path_list.append(abifile_path.split("/")[-1])

        template_aligned, query_aligned_list, nongap, region_start, region_end = recursive_alignment(template, query_list)
         
        if start is None and end is None:
            start = region_start
            end = region_end
        else:
            if region_start > start: 
                nstart = 0 
            else:
                nstart = start - region_start

            if region_end < end:
                nend = region_end - region_start
            else:
                nend = end - region_start 

            m = 0
            n_list = [] 
            for n, char in enumerate(template_aligned.seq):
                if char == "-":
                    pass
                else:
                    if nstart <= m <= nend:
                        n_list.append(n)      
                    else:
                        pass 
                    m += 1
            
            new_abidata_list = []
            new_query_aligned_list = []
            for (aquery, strand), aabidata in zip(query_aligned_list, abidata_list):
                
                m = 0 
                for n, char in enumerate(aquery):
                    if n == min(n_list): 
                        ms = m
                    if n == max(n_list):
                        me = m
                    
                    if char == "-":
                        pass 
                    else:
                        m += 1
                
                new_query_aligned_list.append((aquery[min(n_list):max(n_list)], strand))
                new_abidata = copy.deepcopy(aabidata) 
                
                if strand == 1:
                    nstart = ms
                    nend   = me
                else:
                    nstart, nend   = (len(aabidata["conf"][0]) - me, len(aabidata["conf"][0]) - ms) 
                

                new_abidata["conf"][0] = new_abidata["conf"][0][nstart:nend]
                new_abidata["conf"][1] = new_abidata["conf"][1][nstart:nend] 
                new_abidata["channel"]["A"][0]  = new_abidata["channel"]["A"][0][nstart:nend]
                new_abidata["channel"]["A"][1]  = new_abidata["channel"]["A"][1][nstart:nend] 
                new_abidata["channel"]["T"][0]  = new_abidata["channel"]["T"][0][nstart:nend]
                new_abidata["channel"]["T"][1]  = new_abidata["channel"]["T"][1][nstart:nend] 
                new_abidata["channel"]["G"][0]  = new_abidata["channel"]["G"][0][nstart:nend]
                new_abidata["channel"]["G"][1]  = new_abidata["channel"]["G"][1][nstart:nend] 
                new_abidata["channel"]["C"][0]  = new_abidata["channel"]["C"][0][nstart:nend]
                new_abidata["channel"]["C"][1]  = new_abidata["channel"]["C"][1][nstart:nend] 
                
                for nucl in "ATGC":
                    positions = [p - aabidata["channel"][nucl][0][nstart] for p in aabidata["_channel"][nucl][0] if aabidata["channel"][nucl][0][nstart] <= p <= aabidata["channel"][nucl][0][nend]]
                    values    = [v for p, v in zip(aabidata["_channel"][nucl][0], aabidata["_channel"][nucl][1]) if aabidata["channel"][nucl][0][nstart] <= p <= aabidata["channel"][nucl][0][nend]]
                    new_abidata["_channel"][nucl][0] = positions
                    new_abidata["_channel"][nucl][1] = values
                new_abidata_list.append(new_abidata)

            template_aligned   = cropdna(template_aligned, min(n_list), max(n_list), quinable=0) 
            query_aligned_list = new_query_aligned_list
            abidata_list = new_abidata_list

        i = 0 
        new_abidata_list = [] 
        for query_aligned, strand in query_aligned_list:
            abidata = reform_abidata(abidata_list[i], query_aligned, strand)
            new_abidata_list.append(abidata) 
            i += 1
        
        query_aligned_list = list(list(zip(*query_aligned_list))[0])
        for i in range(len(query_aligned_list)):
            aquery   = query_aligned_list[i]
            
            newquery = ""
            for j, char in enumerate(aquery):
                if char == "-":
                    pass
                else:
                    break
            
            for k, char in enumerate(aquery[::-1]):
                if char == "-":
                    pass
                else:
                    break
            
            if k == 0:
                newquery = j*"/" + aquery[j:]
            else:
                newquery = j*"/" + aquery[j:-1*k] + k*"/"
            
            query_aligned_list[i] = newquery

        if linebreak is None:
            ax_all = visualize(template_aligned, new_abidata_list, query_aligned_list, abiname=abifile_path_list)
        else:
            ax_alls = []
            for i in range(0, len(template_aligned.seq), linebreak):
                ax_all = visualize(template_aligned, new_abidata_list, query_aligned_list, abiname=abifile_path_list, start=i, end=i+linebreak if i+linebreak < len(template_aligned.seq) else len(template_aligned.seq),
                                   display_quality=display_quality)
                if i+linebreak > len(template_aligned.seq):
                    space = (i+linebreak - len(template_aligned.seq)) / linebreak
                    x0, x1, y0, y1 = ax_all.get_inner_corner() 
                    spacer = pw.Brick(ax=pw.basefigure.add_axes([x0, y0, (x1-x0)*(space/(1-space)), y1-y0]))
                    spacer.set_xticks([])
                    spacer.set_yticks([])
                    spacer.spines["right"].set_visible(False)
                    spacer.spines["bottom"].set_visible(False)
                    spacer.spines["left"].set_visible(False)
                    spacer.spines["top"].set_visible(False)
                    spacer.patch.set_alpha(0.0)
                    pw.param["margin"] = None
                    ax_all = ax_all|spacer
                ax_alls.append(ax_all) 
            pw.param["margin"] = 0.4
            ax_all = pw.stack(ax_alls, operator="/")
    else:   
        abidata  = abi_to_dict(abipath)
        query    = generate_consensusseq(abidata)
        result_f = gl_alignment(template, query[0], single=True)
        result_r = gl_alignment(template, query[1], single=True)
        
        if result_f[-1] >= result_r[-1]:
            template_aligned, query_aligned, nongap, region_start, region_end, score = result_f
            strand = 1
        else:
            template_aligned, query_aligned, nongap, region_start, region_end, score = result_r
            strand = -1
        
        if start is None and end is None:
            start = region_start
            end   = region_end
        else:
            if region_start > start: 
                nstart = 0 
            else:
                nstart = start - region_start

            if region_end < end:
                nend = region_end - region_start
            else:
                nend = end - region_start 

            m = 0
            n_list = [] 
            for n, char in enumerate(template_aligned.seq):
                if char == "-":
                    pass
                else:
                    if nstart <= m <= nend:
                        n_list.append(n)      
                    else:
                        pass 
                    m += 1
            
            m = 0 
            for n, char in enumerate(query_aligned):
                if n == min(n_list): 
                    ms = m
                if n == max(n_list):
                    me = m
                
                if char == "-":
                    pass 
                else:
                    m += 1
            
            new_abidata = copy.deepcopy(abidata)  
            if strand == 1:
                nstart = ms
                nend   = me
            else:
                nstart, nend   = (len(abidata["conf"][0]) - me, len(abidata["conf"][0]) - ms) 
            
            new_abidata["conf"][0] = new_abidata["conf"][0][nstart:nend]
            new_abidata["conf"][1] = new_abidata["conf"][1][nstart:nend] 
            new_abidata["channel"]["A"][0]  = new_abidata["channel"]["A"][0][nstart:nend]
            new_abidata["channel"]["A"][1]  = new_abidata["channel"]["A"][1][nstart:nend] 
            new_abidata["channel"]["T"][0]  = new_abidata["channel"]["T"][0][nstart:nend]
            new_abidata["channel"]["T"][1]  = new_abidata["channel"]["T"][1][nstart:nend] 
            new_abidata["channel"]["G"][0]  = new_abidata["channel"]["G"][0][nstart:nend]
            new_abidata["channel"]["G"][1]  = new_abidata["channel"]["G"][1][nstart:nend] 
            new_abidata["channel"]["C"][0]  = new_abidata["channel"]["C"][0][nstart:nend]
            new_abidata["channel"]["C"][1]  = new_abidata["channel"]["C"][1][nstart:nend] 
            
            for nucl in "ATGC":
                positions = [p - abidata["channel"][nucl][0][nstart] for p in abidata["_channel"][nucl][0] if abidata["channel"][nucl][0][nstart] <= p <= abidata["channel"][nucl][0][nend]]
                values    = [v for p, v in zip(abidata["_channel"][nucl][0], abidata["_channel"][nucl][1]) if abidata["channel"][nucl][0][nstart] <= p <= abidata["channel"][nucl][0][nend]]
                new_abidata["_channel"][nucl][0] = positions
                new_abidata["_channel"][nucl][1] = values
            
            template_aligned = cropdna(template_aligned, min(n_list), max(n_list), quinable=0) 
            query_aligned    = query_aligned[min(n_list):max(n_list)]
            abidata          = new_abidata 

        abidata = reform_abidata(abidata, query_aligned, strand)  
        if linebreak is None:
            ax_all  = visualize(template_aligned, abidata, query_aligned, display_quality=display_quality)
        else:
            ax_alls = [] 
            for i in range(0, end-start, linebreak):
                ax_all = visualize(template_aligned, abidata, query_aligned, abiname=abipath.split("/")[-1], start=i, end=i+linebreak if i+linebreak < len(template_aligned.seq) else len(template_aligned.seq), 
                                   display_quality=display_quality)
                if i+linebreak > len(template_aligned.seq):
                    space = (i+linebreak - len(template_aligned.seq)) / linebreak
                    x0, x1, y0, y1 = ax_all.get_inner_corner() 
                    spacer = pw.Brick(ax=pw.basefigure.add_axes([x0, y0, (x1-x0)*(space/(1-space)), y1-y0]))
                    spacer.set_xticks([])
                    spacer.set_yticks([])
                    spacer.spines["right"].set_visible(False)
                    spacer.spines["bottom"].set_visible(False)
                    spacer.spines["left"].set_visible(False)
                    spacer.spines["top"].set_visible(False)
                    spacer.patch.set_alpha(0.0)
                    pw.param["margin"] = None
                    ax_all = ax_all|spacer
                ax_alls.append(ax_all) 
            pw.param["margin"] = 0.4
            ax_all = pw.stack(ax_alls, operator="/")
        
    ax_all.set_suptitle(project + ":" + "{}..{}".format(start, end))
    if output is None:
        pass
    else:
        if dpi is None or output.split(".")[-1] == "pdf":
            ax_all.savefig(output)
        else:
            ax_all.savefig(output, dpi=dpi)
    
    return ax_all 

#if __name__ == "__main__":
#    p = argparse.ArgumentParser()
#    p.add_argument("-q", "--query",     type=str,  help="abi file path or path to the directory containing ab1 files")
#    p.add_argument("-s", "--subject",   type=str,  help="genbank file path")
#    p.add_argument("-l", "--linebreak", type=int,  default=None, help="Sequence length for line break")
#    p.add_argument("-o", "--output",    type=str,  default=None, help="Output file path")
#    p.add_argument("-rs", "--start",    type=int,  default=None, help="Start position of the subject sequence region to be visualized, The output image format should be specified by filename extension.")
#    p.add_argument("-re", "--end",      type=int,  default=None, help="End position of the subject sequence region to be visualized")
#    p.add_argument("-wq", "--quality",  choices=("True", "False"), default="True", help="If True, display bar plot representing Quality value at each nucleotide position.")
#    p.add_argument("-d",  "--dpi",       type=int,  default=None, help="Resolution of the output image. If output format is pdf, the value is ignored.")
#
#    args      = p.parse_args() 
#    abipath   = args.query 
#    gbkpath   = args.subject
#    linebreak = args.linebreak
#    output    = args.output
#    start     = args.start
#    end       = args.end
#    quality   = args.quality
#    quality = True if quality == "True" else False
#    dpi       = args.dpi
#    ax_all = view_sanger(gbkpath, abipath, start, end, linebreak=linebreak, output=output, display_quality=quality, dpi=dpi) 
