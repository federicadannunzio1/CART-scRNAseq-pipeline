import { useState } from "react";

const palette = {
  bg: "#0D1117",
  panel: "#161B22",
  panelBorder: "#21262D",
  accent1: "#E63946",   // CAR scREP gold standard
  accent2: "#2A9D8F",   // kNN / method B
  accent3: "#F4A261",   // method A / signature
  accent4: "#264653",   // both
  green: "#3FB950",
  muted: "#8B949E",
  text: "#E6EDF3",
  subtext: "#8B949E",
  highlight: "#58A6FF",
};

const SUMMARY = [
  { name:"Bo_bone_I",    patient:"Bo", tissue:"Midollo", tp:"I",  n:7300,  tcell:7283, screp:142,  pctScrep:1.95,  nA:429,  sensA:50.7, newA:357, nB:204,  sensB:61.3, newB:117, AB:42  },
  { name:"Ca_bone_I",    patient:"Ca", tissue:"Midollo", tp:"I",  n:1125,  tcell:1125, screp:215,  pctScrep:19.11, nA:78,   sensA:14.9, newA:46,  nB:268,  sensB:35.8, newB:191, AB:15  },
  { name:"Me_bone_I",    patient:"Me", tissue:"Midollo", tp:"I",  n:786,   tcell:786,  screp:78,   pctScrep:9.92,  nA:45,   sensA:11.5, newA:36,  nB:18,   sensB:0,    newB:18,  AB:3   },
  { name:"Ca_blood_AB",  patient:"Ca", tissue:"Sangue",  tp:"AB", n:5353,  tcell:2463, screp:12,   pctScrep:0.22,  nA:62,   sensA:100,  newA:50,  nB:0,    sensB:0,    newB:0,   AB:0   },
  { name:"Ca_bone_AB",   patient:"Ca", tissue:"Midollo", tp:"AB", n:3967,  tcell:null, screp:0,    pctScrep:0,     nA:null, sensA:null, newA:null,nB:null, sensB:null, newB:null,AB:null, skipped:true },
  { name:"Bo_blood_AB",  patient:"Bo", tissue:"Sangue",  tp:"AB", n:6288,  tcell:4610, screp:1566, pctScrep:24.9,  nA:924,  sensA:49.3, newA:153, nB:2138, sensB:92.9, newB:684, AB:124 },
  { name:"Bo_bone_AB",   patient:"Bo", tissue:"Midollo", tp:"AB", n:4275,  tcell:3003, screp:958,  pctScrep:22.41, nA:599,  sensA:52,   newA:103, nB:1285, sensB:91.7, newB:410, AB:73  },
  { name:"Me_bone_AB",   patient:"Me", tissue:"Midollo", tp:"AB", n:5285,  tcell:4601, screp:260,  pctScrep:4.92,  nA:432,  sensA:82.7, newA:217, nB:294,  sensB:81.2, newB:83,  AB:74  },
];

const PATIENT_COLORS = { Bo:"#58A6FF", Ca:"#F4A261", Me:"#3FB950" };
const TP_COLORS = { I:"#6E40C9", AB:"#E63946" };

function SensBar({ val, max=100, color }) {
  if (val === null || val === undefined) return (
    <span style={{color:palette.muted,fontSize:11}}>—</span>
  );
  const w = (val / max) * 80;
  return (
    <div style={{display:"flex",alignItems:"center",gap:6}}>
      <div style={{width:80,height:8,background:"#21262D",borderRadius:4,overflow:"hidden"}}>
        <div style={{width:`${w}%`,height:"100%",background:color,borderRadius:4,
          transition:"width 0.6s ease"}}/>
      </div>
      <span style={{fontSize:11,color:val>70?color:palette.muted,fontWeight:val>70?"700":"400",
        fontFamily:"monospace",minWidth:38}}>{val}%</span>
    </div>
  );
}

function PipelineStep({ num, title, subtitle, color, icon, detail }) {
  return (
    <div style={{display:"flex",flexDirection:"column",alignItems:"center",flex:1,minWidth:0}}>
      <div style={{
        width:56,height:56,borderRadius:"50%",
        background:`${color}22`,border:`2px solid ${color}`,
        display:"flex",alignItems:"center",justifyContent:"center",
        fontSize:22,marginBottom:8,position:"relative"
      }}>
        {icon}
        <div style={{
          position:"absolute",top:-6,right:-6,
          width:20,height:20,borderRadius:"50%",
          background:color,display:"flex",alignItems:"center",
          justifyContent:"center",fontSize:11,fontWeight:700,color:"#0D1117"
        }}>{num}</div>
      </div>
      <div style={{fontWeight:700,color:palette.text,fontSize:13,textAlign:"center",marginBottom:3}}>{title}</div>
      <div style={{fontSize:11,color:palette.subtext,textAlign:"center",lineHeight:1.4}}>{subtitle}</div>
      {detail && (
        <div style={{
          marginTop:8,padding:"4px 8px",background:`${color}11`,
          borderRadius:6,border:`1px solid ${color}44`,
          fontSize:10,color:color,textAlign:"center",fontFamily:"monospace"
        }}>{detail}</div>
      )}
    </div>
  );
}

function Arrow() {
  return (
    <div style={{display:"flex",alignItems:"center",padding:"0 4px",paddingTop:12,flexShrink:0}}>
      <div style={{color:palette.muted,fontSize:20}}>→</div>
    </div>
  );
}

function SensLabel({ v, color }) {
  if (v === null || v === undefined) return <span style={{color:palette.muted}}>—</span>;
  const badge = v >= 70 ? "✓ Buona" : v >= 40 ? "~ Media" : "✗ Bassa";
  const bcolor = v >= 70 ? palette.green : v >= 40 ? "#F4A261" : palette.accent1;
  return (
    <span style={{color:bcolor,fontSize:10,fontWeight:600,
      background:`${bcolor}18`,padding:"1px 5px",borderRadius:3}}>
      {badge} ({v}%)
    </span>
  );
}

export default function App() {
  const [activeRow, setActiveRow] = useState(null);

  const nonSkipped = SUMMARY.filter(d => !d.skipped);

  // Biological highlight data
  const boExpansion = {
    blood_I: 0,
    bone_I: 142,
    blood_AB: 1566,
    bone_AB: 958,
  };

  return (
    <div style={{
      background: palette.bg,
      minHeight: "100vh",
      fontFamily: "'Helvetica Neue', Helvetica, Arial, sans-serif",
      color: palette.text,
      padding: "28px 20px",
      maxWidth: 1100,
      margin: "0 auto",
    }}>

      {/* Header */}
      <div style={{marginBottom:32,borderBottom:`1px solid ${palette.panelBorder}`,paddingBottom:20}}>
        <div style={{display:"flex",alignItems:"flex-start",gap:16}}>
          <div style={{
            width:44,height:44,borderRadius:10,
            background:`${palette.accent1}22`,border:`2px solid ${palette.accent1}`,
            display:"flex",alignItems:"center",justifyContent:"center",
            fontSize:22,flexShrink:0
          }}>🧬</div>
          <div>
            <h1 style={{margin:0,fontSize:22,fontWeight:800,letterSpacing:-0.5,
              background:`linear-gradient(135deg, ${palette.text} 0%, ${palette.highlight} 100%)`,
              WebkitBackgroundClip:"text",WebkitTextFillColor:"transparent"}}>
              CAR-T Cell Detection — Risultati Analisi
            </h1>
            <p style={{margin:"4px 0 0",color:palette.subtext,fontSize:13}}>
              Dataset: 8 campioni • 3 pazienti (Bo, Ca, Me) • 2 timepoint (I, AB) • 2 tessuti (midollo, sangue)
            </p>
          </div>
        </div>
      </div>

      {/* PIPELINE SCHEMA */}
      <div style={{
        background:palette.panel,border:`1px solid ${palette.panelBorder}`,
        borderRadius:14,padding:24,marginBottom:24
      }}>
        <div style={{fontSize:11,fontWeight:700,color:palette.highlight,
          letterSpacing:2,textTransform:"uppercase",marginBottom:16}}>
          Schema della pipeline
        </div>

        <div style={{display:"flex",alignItems:"flex-start",gap:4,flexWrap:"wrap"}}>
          <PipelineStep num="1" icon="📦" color="#58A6FF"
            title="Input"
            subtitle="Oggetto Seurat già annotato (lista 8 campioni)"
            detail="IS_CAR_ALLIN_scREP"
          />
          <Arrow/>
          <PipelineStep num="2" icon="🔬" color="#8957E5"
            title="Subset T cells"
            subtitle="Filtra per cell_type: solo linfociti T (esclude monociti, NK, B...)"
            detail="cell_type ∈ T_CELL_TYPES"
          />
          <Arrow/>
          <PipelineStep num="3" icon="📊" color={palette.accent3}
            title="Metodo A — Firma"
            subtitle="DEG CAR+ vs CAR- → top 50 geni → AddModuleScore → soglia 95° pct"
            detail="wilcox | FDR<0.05 | logFC>0.4"
          />
          <Arrow/>
          <PipelineStep num="4" icon="🕸️" color={palette.accent2}
            title="Metodo B — kNN"
            subtitle="Per ogni T cell: % di k=20 vicini in PCA che sono IS_CAR+ ≥ 0.30"
            detail="PCA 30 dim | k=20"
          />
          <Arrow/>
          <PipelineStep num="5" icon="🎯" color={palette.accent1}
            title="Integrazione A∩B"
            subtitle="5 classi: scREP gold · A+B alta conf. · solo A · solo B · negativo"
            detail="A∩B = massima confidenza"
          />
        </div>

        {/* Legend integration classes */}
        <div style={{
          marginTop:20,padding:"12px 16px",
          background:"#0D111788",borderRadius:10,
          display:"flex",flexWrap:"wrap",gap:12
        }}>
          {[
            {color:"#264653", label:"scREP confirmed", desc:"Gold standard VDJ"},
            {color:palette.accent1, label:"new A∩B", desc:"Alta confidenza"},
            {color:palette.accent3, label:"new A only", desc:"Solo firma — verifica DEG"},
            {color:palette.accent2, label:"new B only", desc:"Solo kNN — verifica UMAP"},
            {color:"#555", label:"CAR negative", desc:"Negativi"},
          ].map(({color,label,desc})=>(
            <div key={label} style={{display:"flex",alignItems:"center",gap:6}}>
              <div style={{width:10,height:10,borderRadius:"50%",background:color,flexShrink:0}}/>
              <span style={{fontSize:11,color:palette.text,fontWeight:600}}>{label}</span>
              <span style={{fontSize:10,color:palette.subtext}}>{desc}</span>
            </div>
          ))}
        </div>
      </div>

      {/* Two boxes: why BCAR3 failed + kNN logic */}
      <div style={{display:"grid",gridTemplateColumns:"1fr 1fr",gap:16,marginBottom:24}}>
        <div style={{
          background:`${palette.accent1}10`,
          border:`1px solid ${palette.accent1}44`,
          borderRadius:12,padding:18
        }}>
          <div style={{fontWeight:700,color:palette.accent1,fontSize:13,marginBottom:8}}>
            ⚠️ Perché l'approccio originale (BCAR3) era sbagliato
          </div>
          <div style={{fontSize:12,color:palette.text,lineHeight:1.7}}>
            Lo script precedente cercava la sottostringa <code style={{background:"#ffffff18",padding:"1px 4px",borderRadius:3,color:palette.accent3}}>CAR</code> nei nomi delle features e selezionava <strong>BCAR3</strong> — un gene endogeno umano (Breast Cancer Anti-Estrogen Resistance 3, chr1), non correlato al costrutto CAR-T. Risultato: overlap sistematicamente 0 su tutti i campioni e sensibilità dello 0%. Non è rumore statistico: le due variabili misuravano fenomeni biologici completamente diversi.
          </div>
        </div>

        <div style={{
          background:`${palette.accent2}10`,
          border:`1px solid ${palette.accent2}44`,
          borderRadius:12,padding:18
        }}>
          <div style={{fontWeight:700,color:palette.accent2,fontSize:13,marginBottom:8}}>
            💡 Perché il Metodo B (kNN) funziona meglio
          </div>
          <div style={{fontSize:12,color:palette.text,lineHeight:1.7}}>
            Le CAR-T si addensano in cluster specifici nello spazio PCA perché condividono un profilo trascrittomica altamente omogeneo (stesso costrutto, stesso differenziamento indotto). Il kNN sfrutta questa struttura geometrica locale: se una cellula è "circondata" da CAR-T confermate nel manifold, è molto probabile che sia anch'essa una CAR-T. Funziona particolarmente bene quando la frequenza delle CAR-T è alta (<strong>Bo_AB: ~23–25%</strong>), formando cluster compatti.
          </div>
        </div>
      </div>

      {/* RESULTS TABLE */}
      <div style={{
        background:palette.panel,border:`1px solid ${palette.panelBorder}`,
        borderRadius:14,padding:24,marginBottom:24
      }}>
        <div style={{fontSize:11,fontWeight:700,color:palette.highlight,
          letterSpacing:2,textTransform:"uppercase",marginBottom:16}}>
          Tabella risultati per campione
        </div>

        <div style={{overflowX:"auto"}}>
          <table style={{width:"100%",borderCollapse:"collapse",fontSize:12}}>
            <thead>
              <tr>
                {["Campione","Paziente","Tessuto","TP","Cellule tot.","T cells","CAR+ scREP","%",
                  "Sens. A","Nuovi A","Sens. B","Nuovi B","A∩B (alta conf.)"
                ].map(h=>(
                  <th key={h} style={{
                    padding:"8px 10px",textAlign:"left",
                    borderBottom:`2px solid ${palette.panelBorder}`,
                    color:palette.subtext,fontWeight:600,
                    fontSize:10,letterSpacing:0.5,whiteSpace:"nowrap"
                  }}>{h}</th>
                ))}
              </tr>
            </thead>
            <tbody>
              {SUMMARY.map((d,i)=>{
                const isActive = activeRow === i;
                const bg = d.skipped ? `${palette.accent1}08`
                  : isActive ? `${palette.highlight}10`
                  : i%2===0 ? "transparent" : `${palette.panelBorder}40`;
                return (
                  <tr key={d.name}
                    onClick={()=>setActiveRow(isActive ? null : i)}
                    style={{background:bg,cursor:"pointer",
                      transition:"background 0.15s",
                      opacity: d.skipped ? 0.5 : 1
                    }}>
                    <td style={{padding:"9px 10px",fontWeight:700,color:palette.text,
                      fontFamily:"monospace",fontSize:11}}>
                      {d.name}
                    </td>
                    <td style={{padding:"9px 10px"}}>
                      <span style={{
                        background:`${PATIENT_COLORS[d.patient]}22`,
                        color:PATIENT_COLORS[d.patient],
                        padding:"2px 8px",borderRadius:20,fontSize:11,fontWeight:700
                      }}>{d.patient}</span>
                    </td>
                    <td style={{padding:"9px 10px",color:palette.subtext,fontSize:11}}>{d.tissue}</td>
                    <td style={{padding:"9px 10px"}}>
                      <span style={{
                        background:`${TP_COLORS[d.tp]}22`,
                        color:TP_COLORS[d.tp],
                        padding:"2px 8px",borderRadius:20,fontSize:11,fontWeight:700
                      }}>{d.tp}</span>
                    </td>
                    <td style={{padding:"9px 10px",color:palette.subtext,fontFamily:"monospace",fontSize:11}}>
                      {d.n.toLocaleString()}
                    </td>
                    <td style={{padding:"9px 10px",color:palette.subtext,fontFamily:"monospace",fontSize:11}}>
                      {d.tcell !== null ? d.tcell.toLocaleString() : "—"}
                    </td>
                    <td style={{padding:"9px 10px"}}>
                      {d.skipped ? (
                        <span style={{color:palette.accent1,fontSize:10,fontWeight:600}}>SKIP (0)</span>
                      ) : (
                        <span style={{
                          color: d.screp > 500 ? palette.accent1 : d.screp > 100 ? palette.accent3 : palette.subtext,
                          fontFamily:"monospace",fontWeight:d.screp>100?"700":"400"
                        }}>{d.screp}</span>
                      )}
                    </td>
                    <td style={{padding:"9px 10px",fontFamily:"monospace",fontSize:11,
                      color: d.pctScrep > 15 ? palette.accent1 : d.pctScrep > 5 ? palette.accent3 : palette.subtext
                    }}>
                      {d.pctScrep > 0 ? `${d.pctScrep}%` : "0%"}
                    </td>
                    <td style={{padding:"9px 10px"}}>
                      <SensBar val={d.sensA} color={palette.accent3}/>
                    </td>
                    <td style={{padding:"9px 10px",fontFamily:"monospace",
                      color:d.newA>200?palette.accent3:palette.subtext,fontSize:11}}>
                      {d.newA !== null ? d.newA : "—"}
                    </td>
                    <td style={{padding:"9px 10px"}}>
                      <SensBar val={d.sensB} color={palette.accent2}/>
                    </td>
                    <td style={{padding:"9px 10px",fontFamily:"monospace",
                      color:d.newB>200?palette.accent2:palette.subtext,fontSize:11}}>
                      {d.newB !== null ? d.newB : "—"}
                    </td>
                    <td style={{padding:"9px 10px"}}>
                      {d.AB !== null ? (
                        <span style={{
                          color: d.AB > 50 ? palette.accent1 : d.AB > 10 ? palette.accent3 : palette.subtext,
                          fontWeight: d.AB > 50 ? "700" : "400",
                          fontFamily:"monospace",fontSize:12
                        }}>{d.AB > 50 ? "⭐ " : ""}{d.AB}</span>
                      ) : <span style={{color:palette.muted}}>—</span>}
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </table>
        </div>
      </div>

      {/* KEY FINDINGS */}
      <div style={{
        background:palette.panel,border:`1px solid ${palette.panelBorder}`,
        borderRadius:14,padding:24,marginBottom:24
      }}>
        <div style={{fontSize:11,fontWeight:700,color:palette.highlight,
          letterSpacing:2,textTransform:"uppercase",marginBottom:20}}>
          Interpretazione biologica e tecnica
        </div>

        <div style={{display:"grid",gridTemplateColumns:"repeat(auto-fit,minmax(300px,1fr))",gap:16}}>

          {/* Finding 1: Bo expansion */}
          <div style={{background:"#0D1117",borderRadius:12,padding:18,
            border:`1px solid ${PATIENT_COLORS.Bo}44`}}>
            <div style={{fontWeight:700,color:PATIENT_COLORS.Bo,fontSize:13,marginBottom:10}}>
              🔷 Paziente Bo — Espansione massiva CAR-T
            </div>
            <div style={{display:"flex",flexDirection:"column",gap:8}}>
              {[
                {label:"Midollo I",val:142,max:1566,color:TP_COLORS.I},
                {label:"Midollo AB",val:958,max:1566,color:TP_COLORS.AB},
                {label:"Sangue AB",val:1566,max:1566,color:TP_COLORS.AB},
              ].map(({label,val,max,color})=>(
                <div key={label} style={{display:"flex",alignItems:"center",gap:8}}>
                  <span style={{fontSize:11,color:palette.subtext,minWidth:80}}>{label}</span>
                  <div style={{flex:1,height:10,background:"#21262D",borderRadius:5,overflow:"hidden"}}>
                    <div style={{
                      width:`${(val/max)*100}%`,height:"100%",
                      background:color,borderRadius:5
                    }}/>
                  </div>
                  <span style={{fontSize:11,fontFamily:"monospace",
                    color:color,fontWeight:700,minWidth:36}}>{val}</span>
                </div>
              ))}
            </div>
            <div style={{marginTop:10,fontSize:11,color:palette.subtext,lineHeight:1.5}}>
              Bo mostra il pattern atteso: incremento netto delle CAR-T tra il timepoint iniziale (I) e il follow-up (AB). Il metodo B funziona molto bene qui (sens. ~92%) perché le CAR-T formano un cluster compatto nel PCA.
            </div>
          </div>

          {/* Finding 2: Sensitivity patterns */}
          <div style={{background:"#0D1117",borderRadius:12,padding:18,
            border:`1px solid ${palette.accent3}44`}}>
            <div style={{fontWeight:700,color:palette.accent3,fontSize:13,marginBottom:10}}>
              📉 Metodo A — Problema di specificità
            </div>
            <div style={{fontSize:11,color:palette.text,lineHeight:1.6}}>
              Il Metodo A tende a <strong style={{color:palette.accent3}}>sovrastimare</strong> i CAR+ in molti campioni:
            </div>
            <div style={{marginTop:10,display:"flex",flexDirection:"column",gap:6}}>
              {[
                {s:"Bo_bone_I",newA:357,screp:142,warn:"357 nuovi vs 142 gold standard"},
                {s:"Ca_bone_I",newA:46,screp:215,warn:"Sensibilità solo 14.9%"},
                {s:"Me_bone_I",newA:36,screp:78,warn:"Sensibilità solo 11.5%"},
              ].map(d=>(
                <div key={d.s} style={{display:"flex",gap:8,alignItems:"flex-start"}}>
                  <span style={{
                    fontFamily:"monospace",fontSize:10,
                    color:palette.accent3,background:`${palette.accent3}18`,
                    padding:"2px 6px",borderRadius:4,flexShrink:0
                  }}>{d.s}</span>
                  <span style={{fontSize:10,color:palette.subtext}}>{d.warn}</span>
                </div>
              ))}
            </div>
            <div style={{marginTop:10,fontSize:11,color:palette.subtext,lineHeight:1.5}}>
              Il DEG probabilmente cattura uno stato di attivazione/effettore generico dei T cells, non specifico del costrutto CAR. Controlla i top DEG: se vedi geni come GZMB, PRF1, IFNG senza marker CAR-specifici è un falso positivo generalizzato.
            </div>
          </div>

          {/* Finding 3: kNN works better */}
          <div style={{background:"#0D1117",borderRadius:12,padding:18,
            border:`1px solid ${palette.accent2}44`}}>
            <div style={{fontWeight:700,color:palette.accent2,fontSize:13,marginBottom:10}}>
              ✅ Metodo B — Performance dipende dalla frequenza
            </div>
            <div style={{display:"flex",flexDirection:"column",gap:6,marginBottom:10}}>
              {nonSkipped.map(d=>(
                <div key={d.name} style={{display:"flex",alignItems:"center",gap:8}}>
                  <span style={{fontSize:10,fontFamily:"monospace",minWidth:90,
                    color:palette.subtext}}>{d.name}</span>
                  <div style={{flex:1,height:7,background:"#21262D",borderRadius:4,overflow:"hidden"}}>
                    <div style={{
                      width:`${d.sensB}%`,height:"100%",
                      background: d.sensB>70?palette.accent2:d.sensB>30?"#F4A261":"#E63946",
                      borderRadius:4
                    }}/>
                  </div>
                  <span style={{fontSize:10,fontFamily:"monospace",minWidth:40,
                    color:d.sensB>70?palette.accent2:palette.subtext}}>{d.sensB}%</span>
                  <span style={{fontSize:9,color:palette.subtext}}>{d.pctScrep}% freq.</span>
                </div>
              ))}
            </div>
            <div style={{fontSize:11,color:palette.subtext,lineHeight:1.5}}>
              Trend chiaro: la sensibilità del kNN scala con la frequenza delle CAR-T nel campione. Sotto ~5% di frequenza (<code style={{color:palette.text}}>Me_bone_I</code>, <code style={{color:palette.text}}>Ca_blood_AB</code>), le cellule CAR+ sono troppo sparse per formare un cluster rilevabile.
            </div>
          </div>

          {/* Finding 4: Ca bone AB */}
          <div style={{background:"#0D1117",borderRadius:12,padding:18,
            border:`1px solid #E6D74944`}}>
            <div style={{fontWeight:700,color:"#E6D749",fontSize:13,marginBottom:10}}>
              ❓ Ca_bone_AB — Anomalia da verificare
            </div>
            <div style={{fontSize:11,color:palette.text,lineHeight:1.6}}>
              Questo campione ha <strong style={{color:palette.accent1}}>0 cellule IS_CAR_ALLIN_scREP = YES</strong>
              , quindi è stato skippato. Possibili cause:
            </div>
            <div style={{marginTop:8,display:"flex",flexDirection:"column",gap:5}}>
              {[
                "File VDJ non processato o non linkato",
                "Recovery VDJ molto bassa in questo campione",
                "Vera assenza di CAR-T nel midollo di Ca al timepoint AB (biologicamente possibile)",
                "Bug nella pipeline scREP upstream",
              ].map((item,i)=>(
                <div key={i} style={{display:"flex",gap:6,fontSize:11,color:palette.subtext}}>
                  <span style={{color:"#E6D749",flexShrink:0}}>•</span>
                  <span>{item}</span>
                </div>
              ))}
            </div>
            <div style={{marginTop:8,padding:"8px 10px",
              background:"#E6D74910",borderRadius:7,
              fontSize:11,color:"#E6D749"}}>
              Confrontare con Ca_blood_AB (12 cellule, 0.22%): Ca ha poche CAR-T nel sangue al timepoint AB — biologicamente coerente con una scarsa espansione rispetto a Bo.
            </div>
          </div>

        </div>
      </div>

      {/* RELIABILITY GUIDE */}
      <div style={{
        background:palette.panel,border:`1px solid ${palette.panelBorder}`,
        borderRadius:14,padding:24,marginBottom:16
      }}>
        <div style={{fontSize:11,fontWeight:700,color:palette.highlight,
          letterSpacing:2,textTransform:"uppercase",marginBottom:16}}>
          Come leggere i risultati — Guida alla confidenza
        </div>

        <div style={{display:"flex",flexDirection:"column",gap:10}}>
          {[
            {
              cat:"⭐ scREP_confirmed",
              color:palette.accent4,
              conf:"Massima",
              desc:"Gold standard VDJ. Prova diretta della presenza del TCR/BCR del costrutto.",
              use:"Usa per tutti gli downstream: espansione clonale, differenziamento, funzione."
            },
            {
              cat:"🔴 new_A_and_B",
              color:palette.accent1,
              conf:"Alta",
              desc:"Positivi a entrambi i metodi indipendenti. La convergenza riduce molto i falsi positivi.",
              use:"Utile come supplemento al gold standard. Verifica su UMAP che siano nel cluster CAR+."
            },
            {
              cat:"🟠 new_A_only",
              color:palette.accent3,
              conf:"Bassa — da verificare",
              desc:"Solo la firma trascrittomica. Alta probabilità di T cells effettori endogeni con fenotipo simile.",
              use:"Controlla i top DEG: se dominati da geni di attivazione generica (GZMB, PRF1) sono probabilmente falsi positivi."
            },
            {
              cat:"🟢 new_B_only",
              color:palette.accent2,
              conf:"Media — dipende dal contesto",
              desc:"Solo vicinanza nel PCA alle CAR-T note. Possibili T cells endogeni nello stesso cluster.",
              use:"Controlla su UMAP: se le cellule sono fisicamente adiacenti ai scREP_confirmed sono candidati affidabili."
            },
          ].map(({cat,color,conf,desc,use})=>(
            <div key={cat} style={{
              display:"flex",gap:14,padding:"12px 14px",
              background:"#0D111788",borderRadius:9,
              borderLeft:`3px solid ${color}`
            }}>
              <div style={{minWidth:120,flexShrink:0}}>
                <div style={{fontWeight:700,color:color,fontSize:12}}>{cat}</div>
                <div style={{fontSize:10,color:palette.subtext,marginTop:3}}>
                  Confidenza: <span style={{color:color}}>{conf}</span>
                </div>
              </div>
              <div style={{flex:1}}>
                <div style={{fontSize:11,color:palette.text,marginBottom:4}}>{desc}</div>
                <div style={{fontSize:10,color:palette.subtext}}>
                  <strong style={{color:palette.muted}}>Uso consigliato:</strong> {use}
                </div>
              </div>
            </div>
          ))}
        </div>
      </div>

      {/* Footer note */}
      <div style={{fontSize:10,color:palette.subtext,textAlign:"center",paddingTop:8,
        borderTop:`1px solid ${palette.panelBorder}`,lineHeight:1.6}}>
        ⚠️ <strong style={{color:palette.muted}}>Limite invalicabile:</strong> senza sequenze del costrutto nei FASTQ o VDJ recovery, nessun metodo trascrittomica può escludere T cells endogeni con fenotipo effettore simile alle CAR-T. La colonna <code style={{color:palette.text}}>new_A_and_B</code> è la stima più conservativa ottenibile da questo dataset.
      </div>

    </div>
  );
}
