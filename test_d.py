import os
import re
import json
import urllib.request
import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import networkx as nx

# Configuração da página
st.set_page_config(page_title="Consolidado por Gene", layout="wide")

URL_BASE_ALPHAFOLD = "https://alphafold.ebi.ac.uk/entry/"

MAPA_TRADUCAO_CLINICA = {
    'uncertain significance': 'Sig. Incerta',
    'likely pathogenic': 'Prov. Patogênica',
    'pathogenic': 'Patogênica',
    'likely benign': 'Prov. Benigna',
    'benign': 'Benigna',
    'conflicting interpretations of pathogenicity': 'Interpr. Conflitantes',
    'not provided': 'Não Informado'
}

# --- ESTILIZAÇÃO CUSTOMIZADA PREMIUM (CSS Leve & Arredondado) ---
st.markdown("""
<style>
    /* Fundo geral e tipografia */
    .stApp {
        background-color: #0B0F17;
        font-family: 'Inter', system-ui, -apple-system, sans-serif;
    }

    header {visibility: hidden;}
    .block-container {
        padding-top: 1.5rem;
        padding-bottom: 2rem;
    }

    /* --- CABEÇALHO HERO MODERNO --- */
    .hero-header {
        background: linear-gradient(135deg, #131C2E 0%, #0F172A 100%);
        border: 1px solid #1E293B;
        border-radius: 20px;
        padding: 24px 30px;
        margin-bottom: 24px;
        box-shadow: 0 10px 30px -10px rgba(0, 0, 0, 0.5);
        display: flex;
        align-items: center;
        justify-content: space-between;
    }
    .hero-title-group {
        display: flex;
        align-items: center;
        gap: 16px;
    }
    .hero-title {
        font-size: 2.4rem;
        font-weight: 800;
        letter-spacing: -0.5px;
        background: linear-gradient(90deg, #60A5FA 0%, #A78BFA 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        margin: 0;
    }
    .hero-subtitle {
        color: #64748B;
        font-size: 0.9rem;
        margin-top: 2px;
    }

    /* CARD DE INFORMAÇÕES DO GENE */
    .header-card {
        background: #111827;
        border: 1px solid #1F2937;
        border-radius: 18px;
        padding: 20px 24px;
        margin-bottom: 20px;
    }
    .gene-title {
        font-size: 1.7rem;
        font-weight: 700;
        color: #F9FAFB;
        margin: 0;
        display: flex;
        align-items: center;
        gap: 12px;
    }
    .uniprot-badge {
        background: rgba(16, 185, 129, 0.12);
        color: #10B981;
        border: 1px solid rgba(16, 185, 129, 0.3);
        padding: 5px 14px;
        border-radius: 9999px;
        font-family: monospace;
        font-size: 0.85rem;
        font-weight: 600;
    }

    /* METRIC CARDS ARREDONDADOS */
    div[data-testid="stMetric"] {
        background: #111827;
        border: 1px solid #1F2937;
        padding: 16px 20px;
        border-radius: 16px;
        transition: transform 0.2s ease, border-color 0.2s ease;
    }
    div[data-testid="stMetric"]:hover {
        transform: translateY(-3px);
        border-color: #374151;
    }
    div[data-testid="stMetricLabel"] {
        color: #9CA3AF !important;
        font-size: 0.78rem !important;
        font-weight: 600 !important;
        text-transform: uppercase;
    }

    /* --- ABAS (TABS) MODERNAS & ARREDONDADAS --- */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
        background-color: #111827;
        padding: 8px;
        border-radius: 16px;
        border: 1px solid #1F2937;
    }
    .stTabs [data-baseweb="tab"] {
        height: 44px;
        border-radius: 10px;
        color: #9CA3AF;
        font-weight: 500;
        border: none !important;
        padding: 0 18px;
        transition: all 0.2s ease;
    }
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, #2563EB 0%, #1D4ED8 100%) !important;
        color: #FFFFFF !important;
        font-weight: 600;
        box-shadow: 0 4px 12px rgba(37, 99, 235, 0.3);
    }

    /* --- DATAFRAME/TABELA ARREDONDADA E MODERNA --- */
    div[data-testid="stDataFrame"] {
        background-color: #111827;
        border: 1px solid #1F2937;
        border-radius: 16px;
        padding: 6px;
        overflow: hidden;
    }
</style>
""", unsafe_allow_html=True)

def traduzir_significancia(val):
    if pd.isna(val):
        return val
    val_str = str(val).strip().lower()
    for eng, pt in MAPA_TRADUCAO_CLINICA.items():
        val_str = val_str.replace(eng, pt)
    val_str = val_str.title()
    if len(val_str) > 25:
        return "Conflitante / Múltipla"
    return val_str

@st.cache_data(show_spinner="Carregando base de dados...")
def carregar_dados_otimizado():
    dir_base = os.path.dirname(os.path.abspath(__file__))
    caminhos = [
        os.path.join(dir_base, "resultados", "dbsnp_alphamissense_unificado.tsv"),
        os.path.join(dir_base, "comparacao_dbsnp_alphamissense.tsv"),
        os.path.join(dir_base, "tabela_comparativa_completa.tsv"),
        os.path.join(dir_base, "dbsnp_missense_resultados.tsv"),
        os.path.join(dir_base, "alphamissense_resultados.tsv"),
        os.path.join(dir_base, "recuperados_uniprot.tsv")
    ]
    caminho_escolhido = next((c for c in caminhos if os.path.exists(c)), None)
    if not caminho_escolhido:
        return pd.DataFrame()

    df = pd.read_csv(caminho_escolhido, sep="\t", on_bad_lines="skip", low_memory=False)
    df.columns = df.columns.str.strip()

    if "Gene_Original" in df.columns and "Gene" not in df.columns:
        df.rename(columns={"Gene_Original": "Gene"}, inplace=True)

    col_rs = next((c for c in df.columns if 'rs' in c.lower() or 'dbsnp' in c.lower()), None)
    col_change = next((c for c in df.columns if any(k in c.lower() for k in ["proteina_mudanca", "protein_change", "aa_change"])), None)
    
    if col_rs and col_change:
        df.drop_duplicates(subset=[col_rs, col_change], keep='first', inplace=True)

    col_class = next((c for c in df.columns if any(k in c.lower() for k in ["classificacao", "clinical", "significancia", "clinvar"])), None)
    if col_class:
        df[col_class] = df[col_class].apply(traduzir_significancia)
        
    return df

@st.cache_data(ttl=86400, show_spinner=False)
def buscar_info_uniprot_api(uniprot_id):
    """ Busca FASTA e anotações funcionais reais via REST UniProt """
    if not uniprot_id or str(uniprot_id).lower() in ["n/a", "none", "nan"]:
        return None, {}
    
    fasta_seq = None
    anotacoes = {"signal_peptide": [], "signal_range": None, "domains": [], "binding_sites": [], "interactions": []}
    
    try:
        url_fasta = f"https://www.uniprot.org/uniprot/{uniprot_id}.fasta"
        req = urllib.request.Request(url_fasta, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=5) as resp:
            fasta_seq = resp.read().decode('utf-8')
    except Exception:
        fasta_seq = None

    try:
        url_json = f"https://www.uniprot.org/uniprot/{uniprot_id}.json"
        req = urllib.request.Request(url_json, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=5) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            features = data.get("features", [])
            for f in features:
                ftype = f.get("type")
                desc = f.get("description", "")
                location = f.get("location", {})
                start = location.get("start", {}).get("value")
                end = location.get("end", {}).get("value")
                pos_str = f"({start}-{end})" if start and end else ""

                if ftype == "Signal":
                    anotacoes["signal_peptide"].append(f"Peptídeo Sinal {pos_str} {desc}".strip())
                    if start and end:
                        anotacoes["signal_range"] = (int(start), int(end))
                elif ftype in ["Domain", "Region", "Repeat"]:
                    if desc:
                        anotacoes["domains"].append(f"{desc} {pos_str}".strip())
                elif ftype in ["Binding site", "Active site", "Site", "Metal binding"]:
                    if desc:
                        anotacoes["binding_sites"].append(f"{ftype} {pos_str}: {desc}".strip())
                    
            comments = data.get("comments", [])
            for c in comments:
                if c.get("commentType") == "INTERACTION":
                    for inter in c.get("interactions", []):
                        partner_gene = inter.get("partner1", {}).get("geneName") or inter.get("partner2", {}).get("geneName")
                        partner_acc = inter.get("partner2", {}).get("uniProtKBAccession")
                        if partner_gene:
                            anotacoes["interactions"].append(f"Gene: {partner_gene}")
                        elif partner_acc:
                            anotacoes["interactions"].append(f"UniProt: {partner_acc}")
    except Exception:
        pass

    return fasta_seq, anotacoes

@st.cache_data(ttl=86400, show_spinner="Baixando estrutura 3D...")
def obter_pdb_alphafold(uniprot_id):
    """ Busca o arquivo PDB usando a API oficial da AlphaFold DB """
    if not uniprot_id or str(uniprot_id).lower() in ["n/a", "none", "nan"]:
        return None
    
    uniprot_clean = str(uniprot_id).strip().upper()

    try:
        url_api = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_clean}"
        req = urllib.request.Request(url_api, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=8) as resp:
            data = json.loads(resp.read().decode('utf-8'))
            if isinstance(data, list) and len(data) > 0:
                pdb_url = data[0].get("pdbUrl")
                if pdb_url:
                    req_pdb = urllib.request.Request(pdb_url, headers={'User-Agent': 'Mozilla/5.0'})
                    with urllib.request.urlopen(req_pdb, timeout=10) as resp_pdb:
                        return resp_pdb.read().decode('utf-8')
    except Exception:
        pass

    urls_fallback = [
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_clean}-F1-model_v4.pdb",
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_clean}-F1-model_v3.pdb",
        f"https://files.rcsb.org/download/{uniprot_clean}.pdb"
    ]

    for url in urls_fallback:
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=8) as resp:
                conteudo = resp.read().decode('utf-8')
                if "ATOM" in conteudo:
                    return conteudo
        except Exception:
            continue

    return None

def extrair_coordenadas_ca(pdb_text):
    """ Extrai resíduos e coordenadas 3D de C-Alpha do arquivo PDB """
    ca_atoms = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_num = int(line[22:26].strip())
            res_name = line[17:20].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            ca_atoms.append({"res_num": res_num, "res_name": res_name, "x": x, "y": y, "z": z})
    return pd.DataFrame(ca_atoms)

def calcular_rede_ligacoes(df_ca, distancia_corte=8.0):
    """ Calcula rede de contatos/ligações entre resíduos baseada na distância espacial """
    if df_ca.empty:
        return nx.Graph(), 0, pd.DataFrame()
    
    coords = df_ca[['x', 'y', 'z']].values
    n = len(df_ca)
    
    G = nx.Graph()
    for _, row in df_ca.iterrows():
        G.add_node(row['res_num'], name=f"{row['res_name']}{row['res_num']}")

    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dist_matrix = np.sqrt(np.sum(diff ** 2, axis=-1))

    for i in range(n):
        for j in range(i + 2, n):
            d = dist_matrix[i, j]
            if d <= distancia_corte:
                r1, r2 = df_ca.iloc[i]['res_num'], df_ca.iloc[j]['res_num']
                G.add_edge(r1, r2, weight=d)

    grados = dict(G.degree())
    df_ligacoes = pd.DataFrame([
        {"Posição (AA)": node, "Resíduo": G.nodes[node]['name'], "Número de Ligações (3D)": deg}
        for node, deg in grados.items()
    ]).sort_values(by="Número de Ligações (3D)", ascending=False)

    return G, G.number_of_edges(), df_ligacoes

def limpar_e_converter(serie):
    s = serie.astype(str).str.replace(',', '.').str.strip()
    s = s.replace(['nan', 'none', 'n/a', '-', 'null', '', '<na>', 'None', 'NONE'], np.nan)
    return pd.to_numeric(s, errors='coerce')

def extrair_posicao_de_mudanca(serie):
    def _extrai_pos(x):
        if not isinstance(x, str) or x.lower() in ['none', 'nan', '', 'null', 'n/a']:
            return np.nan
        m = re.search(r"\d+", x)
        return float(m.group(0)) if m else np.nan
    return serie.apply(_extrai_pos)

def extrair_detalhes_proteina(serie):
    def _extrai(x):
        if not isinstance(x, str) or x.lower() in ['none', 'nan', '', 'null', 'n/a']:
            return pd.Series([np.nan, np.nan, np.nan])
        
        m = re.search(r"(?:p\.)?([A-Za-z]{1,3})?(\d+)([A-Za-z]{1,3})?", x)
        if m:
            aa_ref = m.group(1) if m.group(1) else np.nan
            pos = float(m.group(2))
            aa_alt = m.group(3) if m.group(3) else np.nan
            return pd.Series([pos, aa_ref, aa_alt])
        return pd.Series([np.nan, np.nan, np.nan])

    res = serie.apply(_extrai)
    res.columns = ['Posicao_Proteina', 'AA_Referencia', 'AA_Alterado']
    return res

# Inicialização
df = carregar_dados_otimizado()

st.title("🧬 Ren")

if df.empty:
    st.error("Nenhum arquivo de dados consolidado foi encontrado.")
else:
    col_gene = next((c for c in df.columns if c.lower() in ['gene', 'gene_name', 'symbol']), None)

    if not col_gene:
        st.error("Coluna de identificação do Gene não foi encontrada nos dados.")
    else:
        genes_disponiveis = sorted(df[col_gene].dropna().astype(str).str.strip().str.upper().unique())
        
        gene_selecionado = st.selectbox(
            "🔎 Busque ou selecione o Gene:",
            options=[""] + genes_disponiveis,
            index=0,
            placeholder="Digite para buscar (Ex: ABCB11, OLR1, TP53...)"
        )

        if gene_selecionado:
            df_gene = df[df[col_gene].astype(str).str.strip().str.upper() == gene_selecionado].copy()

            if not df_gene.empty:
                col_uniprot = next((c for c in df_gene.columns if any(k in c.lower() for k in ['uniprot_id', 'id_uniprot', 'uniprot', 'accession'])), None)
                col_nome_prot = next((c for c in df_gene.columns if any(k in c.lower() for k in ['nome_proteina', 'protein_name', 'proteina'])), None)
                col_rs = next((c for c in df_gene.columns if 'rs' in c.lower() or 'dbsnp' in c.lower()), None)
                col_score = next((c for c in df_gene.columns if any(k in c.lower() for k in ['alphamissense_score', 'am_pathogenicity', 'score', 'pathogenicity'])), None)
                col_class_dbsnp = next((c for c in df_gene.columns if any(k in c.lower() for k in ["classificacao", "clinical", "significancia", "clinvar"])), None)
                col_pos = next((c for c in df_gene.columns if any(k in c.lower() for k in ['posicao_proteina', 'pos_proteina', 'protein_pos', 'aa_pos'])), None)
                col_change = next((c for c in df_gene.columns if any(k in c.lower() for k in ["proteina_mudanca", "protein_change", "aa_change"])), None)

                uniprot_val = "N/A"
                if col_uniprot:
                    valid_uniprots = df_gene[col_uniprot].dropna().astype(str).str.strip()
                    valid_uniprots = valid_uniprots[~valid_uniprots.str.lower().isin(['none', 'nan', '', 'null', 'n/a'])]
                    if not valid_uniprots.empty:
                        uniprot_val = valid_uniprots.iloc[0]

                nome_prot_val = "N/A"
                if col_nome_prot:
                    valid_prots = df_gene[col_nome_prot].dropna().astype(str).str.strip()
                    valid_prots = valid_prots[~valid_prots.str.lower().isin(['none', 'nan', '', 'null', 'n/a'])]
                    if not valid_prots.empty:
                        nome_prot_val = valid_prots.iloc[0]

                if col_change:
                    detalhes_aa = extrair_detalhes_proteina(df_gene[col_change])
                    df_gene['Posicao_Proteina'] = detalhes_aa['Posicao_Proteina']
                    df_gene['AA_Ref'] = detalhes_aa['AA_Referencia']
                    df_gene['AA_Alt'] = detalhes_aa['AA_Alterado']
                    df_gene['Mudanca_Proteina'] = df_gene[col_change]
                else:
                    df_gene['Posicao_Proteina'] = np.nan
                    df_gene['AA_Ref'] = np.nan
                    df_gene['AA_Alt'] = np.nan
                    df_gene['Mudanca_Proteina'] = np.nan

                col_gene_pos = next((c for c in df_gene.columns if any(k in c.lower() for k in ['posicao_gene', 'gene_pos', 'cds_pos', 'pos_gene', 'hgvs_c', 'cdna'])), None)
                col_genomic_pos = next((c for c in df_gene.columns if any(k in c.lower() for k in ['posicao_genomica', 'genomic_pos', 'chr_pos', 'hgvs_g', 'chrom'])), None)
                col_consequence = next((c for c in df_gene.columns if any(k in c.lower() for k in ['consequencia', 'consequence', 'so_term', 'functional_class'])), None)

                df_gene['Posicao_Gene'] = df_gene[col_gene_pos] if col_gene_pos else np.nan
                df_gene['Posicao_Genomica'] = df_gene[col_genomic_pos] if col_genomic_pos else np.nan
                df_gene['Consequencia_Gene'] = df_gene[col_consequence] if col_consequence else np.nan
                df_gene['Link_AlphaFold'] = f"{URL_BASE_ALPHAFOLD}{uniprot_val}" if uniprot_val != "N/A" else None
                df_gene['AM_Score_Num'] = limpar_e_converter(df_gene[col_score]) if col_score else np.nan

                fasta_original, anotacoes_uniprot = buscar_info_uniprot_api(uniprot_val)

                # --- HEADER CARD ---
                st.markdown(f"""
                <div class="header-card">
                    <div class="gene-title">
                        📌 Gene: <span style="color: #60A5FA;">{gene_selecionado}</span>
                        <span class="uniprot-badge">UniProt: {uniprot_val}</span>
                    </div>
                    {"<div class='protein-subtitle'><b>Proteína:</b> " + nome_prot_val + "</div>" if nome_prot_val != "N/A" else ""}
                </div>
                """, unsafe_allow_html=True)

                # --- METRIC CARDS ---
                c1, c2, c3, c4 = st.columns(4)
                with c1:
                    st.metric("Total de Registros", f"{len(df_gene):,}".replace(",", "."))
                with c2:
                    if col_rs:
                        rs_vals = df_gene[col_rs].dropna().astype(str).str.strip()
                        rs_count = rs_vals[~rs_vals.isin(['-', 'nan', '', 'N/A', 'None', '<NA>'])].nunique()
                        st.metric("Variantes dbSNP (rsID)", f"{rs_count:,}".replace(",", "."))
                    else:
                        st.metric("Variantes dbSNP", "0")
                with c3:
                    trocas = df_gene['Mudanca_Proteina'].dropna().astype(str).str.strip()
                    trocas_validas = trocas[~trocas.isin(['-', 'nan', '', 'N/A', 'None', '<NA>', 'null'])].nunique()
                    st.metric("Trocas Únicas (AA)", f"{trocas_validas:,}".replace(",", "."))
                with c4:
                    vals_s = df_gene['AM_Score_Num'].dropna()
                    m_score = vals_s.mean() if not vals_s.empty else np.nan
                    st.metric("Score Médio AlphaMissense", f"{m_score:.3f}" if not np.isnan(m_score) else "N/A")

                st.divider()

                # PAINEL PLOTLY
                st.subheader("📊 Painel Visual da Proteína")
                df_valid = df_gene.dropna(subset=['AM_Score_Num']).copy()

                if not df_valid.empty:
                    fig = make_subplots(
                        rows=2, cols=2,
                        specs=[[{"colspan": 2}, None], [{}, {}]],
                        subplot_titles=(
                            "Pontuação de Patogenocidade por Posição na Proteína",
                            "Distribuição Numérica dos Scores",
                            "Validação: Scores vs Classificação dbSNP / ClinVar"
                        ),
                        vertical_spacing=0.18,
                        horizontal_spacing=0.1
                    )

                    eixo_x = df_valid['Posicao_Proteina'] if df_valid['Posicao_Proteina'].notna().any() else np.arange(len(df_valid))
                    
                    fig.add_trace(
                        go.Scatter(
                            x=eixo_x,
                            y=df_valid['AM_Score_Num'],
                            mode='markers',
                            marker=dict(size=7, color='#2b5c8f', opacity=0.7),
                            name='Variantes',
                            text=df_valid[col_rs] if col_rs else None,
                            customdata=df_valid['Mudanca_Proteina'] if 'Mudanca_Proteina' in df_valid.columns else None,
                            hovertemplate="<b>rsID:</b> %{text}<br><b>Mudança:</b> %{customdata}<br><b>Posição AA:</b> %{x}<br><b>Score:</b> %{y:.3f}<extra></extra>"
                        ),
                        row=1, col=1
                    )

                    fig.add_hline(y=0.56, line_dash="dash", line_color="red", annotation_text="Patogênico (>= 0.56)", annotation_position="top left", row=1, col=1)
                    fig.add_hline(y=0.34, line_dash="dash", line_color="green", annotation_text="Benigno (<= 0.34)", annotation_position="bottom left", row=1, col=1)

                    fig.add_trace(
                        go.Histogram(x=df_valid['AM_Score_Num'], nbinsx=15, marker_color='#4682b4', name='Frequência'),
                        row=2, col=1
                    )

                    if col_class_dbsnp and col_class_dbsnp in df_valid.columns:
                        df_box = df_valid.dropna(subset=[col_class_dbsnp])
                        for cat in df_box[col_class_dbsnp].unique():
                            sub_df = df_box[df_box[col_class_dbsnp] == cat]
                            fig.add_trace(
                                go.Box(y=sub_df['AM_Score_Num'], name=str(cat), boxpoints='all', jitter=0.3, pointpos=-1.8, marker=dict(size=4)),
                                row=2, col=2
                            )

                    fig.update_layout(height=600, showlegend=False, template="plotly_white", margin=dict(l=20, r=20, t=50, b=40))
                    fig.update_xaxes(title_text="Posição do Aminoácido", row=1, col=1)
                    fig.update_yaxes(title_text="AlphaMissense Score", row=1, col=1)
                    fig.update_xaxes(title_text="AlphaMissense Score", row=2, col=1)
                    fig.update_yaxes(title_text="Frequência", row=2, col=1)
                    fig.update_xaxes(title_text="Classificação", tickangle=-25, tickfont=dict(size=10), row=2, col=2)
                    fig.update_yaxes(title_text="AlphaMissense Score", row=2, col=2)

                    event = st.plotly_chart(fig, use_container_width=True, on_select="rerun")

                    if event and "select" in event and event["select"]["points"]:
                        pt = event["select"]["points"][0]
                        pos_clicada = pt.get("x")
                        
                        if pos_clicada is not None:
                            df_detalhe = df_valid[df_valid['Posicao_Proteina'] == pos_clicada]

                            if not df_detalhe.empty:
                                st.info(f"🎯 **Sítio Selecionado no Gráfico: Posição `{int(pos_clicada)}`**")
                                cols_det = [col_rs, 'Mudanca_Proteina', 'Posicao_Proteina', col_score, col_class_dbsnp]
                                cols_existentes = [c for c in cols_det if c in df_detalhe.columns]
                                
                                st.dataframe(
                                    df_detalhe[cols_existentes],
                                    use_container_width=True,
                                    hide_index=True,
                                    column_config={
                                        col_score: st.column_config.NumberColumn("AlphaMissense Score", format="%.3f"),
                                        "Posicao_Proteina": st.column_config.NumberColumn("Posição (AA)", format="%d")
                                    }
                                )

                st.divider()

                # MUTAÇÕES POR SÍTIO
                if col_pos and df_gene[col_pos].notna().any():
                    posicoes = limpar_e_converter(df_gene[col_pos])
                elif 'Posicao_Proteina' in df_gene.columns and df_gene['Posicao_Proteina'].notna().any():
                    posicoes = df_gene['Posicao_Proteina']
                elif col_change and col_change in df_gene.columns:
                    posicoes = extrair_posicao_de_mudanca(df_gene[col_change])
                else:
                    posicoes = pd.Series(dtype=float)

                posicoes_validas = posicoes.dropna()
                st.subheader(f"🧬 Mutações por Sítio ({len(posicoes_validas)})")

                if not posicoes_validas.empty:
                    contagem_sitio = posicoes_validas.astype(int).value_counts().sort_index()

                    fig_sitio = go.Figure()

                    for pos, qtd in contagem_sitio.items():
                        cor_haste = '#FF4B5C' if qtd > 1 else '#4A5568'
                        fig_sitio.add_trace(go.Scatter(
                            x=[pos, pos],
                            y=[0, qtd],
                            mode='lines',
                            line=dict(color=cor_haste, width=1.5 if qtd == 1 else 2.5),
                            hoverinfo='skip',
                            showlegend=False
                        ))

                    fig_sitio.add_trace(go.Scatter(
                        x=contagem_sitio.index,
                        y=contagem_sitio.values,
                        mode='markers',
                        marker=dict(
                            size=[10 if y > 1 else 5 for y in contagem_sitio.values],
                            color=['#FF4B5C' if y > 1 else '#718096' for y in contagem_sitio.values],
                            line=dict(width=1, color='#FFFFFF')
                        ),
                        hovertemplate="<b>Sítio (Posição):</b> %{x}<br><b>Número de Mutações:</b> %{y}<extra></extra>",
                        showlegend=False
                    ))

                    fig_sitio.update_layout(
                        paper_bgcolor='rgba(0,0,0,0)',
                        plot_bgcolor='rgba(0,0,0,0)',
                        height=400,
                        margin=dict(l=20, r=20, t=30, b=20),
                        xaxis=dict(
                            title="Posição na proteína (sítio)",
                            showgrid=False,
                            zeroline=True,
                            zerolinecolor='#4A5568'
                        ),
                        yaxis=dict(
                            title="Número de mutações",
                            dtick=1,
                            gridcolor='rgba(255, 255, 255, 0.1)',
                            zeroline=True,
                            zerolinecolor='#4A5568'
                        )
                    )

                    st.caption(f"{len(posicoes_validas)} variantes distribuídas em {contagem_sitio.shape[0]} sítios diferentes.")
                    st.plotly_chart(fig_sitio, use_container_width=True)
                    
                else:
                    st.info("Nenhuma posição válida encontrada para este gene.")

                st.divider()

                # ESTRUTURA POR ABAS
                st.subheader("📑 Detalhamento das Alterações")
                tab_proteina, tab_gene, tab_resumo, tab_3d = st.tabs([
                    "🧬 Alterações na Proteína",
                    "🧬 Alterações no Gene, Genômica & FASTA",
                    "📋 Resumo dos Atributos & Anotações",
                    "🧊 Estrutura 3D, Sinalização de Peptídeo Sinal e Rede de Contatos Espaciais"
                ])

                colunas_redundantes = ['AlphaMissense_Classificacao', 'AM_Score_Num']

                with tab_proteina:
                    st.markdown("##### 🔍 Mutações com Foco na Estrutura Proteica e Aminoácidos")
                    cols_preferenciais = [col_rs, 'Mudanca_Proteina', 'Posicao_Proteina', 'AA_Ref', 'AA_Alt', col_score, col_class_dbsnp, 'Link_AlphaFold']
                    cols_finais_prot = [c for c in cols_preferenciais if c and c in df_gene.columns]
                    outras = [c for c in df_gene.columns if c not in cols_finais_prot and c not in colunas_redundantes + ['Posicao_Gene', 'Posicao_Genomica', 'Consequencia_Gene']]

                    st.dataframe(
                        df_gene[cols_finais_prot + outras],
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            col_score: st.column_config.NumberColumn("AlphaMissense Score", format="%.3f"),
                            col_class_dbsnp: st.column_config.TextColumn("Significância Clínica"),
                            "Mudanca_Proteina": st.column_config.TextColumn("Mudança na Proteína (HGVS p.)"),
                            "Posicao_Proteina": st.column_config.NumberColumn("Posição (AA)", format="%d"),
                            "AA_Ref": st.column_config.TextColumn("AA Ref"),
                            "AA_Alt": st.column_config.TextColumn("AA Alt"),
                            "Link_AlphaFold": st.column_config.LinkColumn("AlphaFold 3D", display_text="Ver Modelo 3D")
                        }
                    )

                with tab_gene:
                    st.markdown("##### 🧬 Mutações com Foco em Posição no Gene (cDNA/CDS), Genômica e Consequência")
                    cols_pref_gene = [col_rs, 'Consequencia_Gene', 'Mudanca_Proteina', col_score, col_class_dbsnp, 'Posicao_Gene', 'Posicao_Genomica']
                    cols_finais_gene = [c for c in cols_pref_gene if c and c in df_gene.columns]
                    outras_gene = [c for c in df_gene.columns if c not in cols_finais_gene and c not in colunas_redundantes + ['Link_AlphaFold', 'Posicao_Proteina', 'AA_Ref', 'AA_Alt']]

                    st.dataframe(
                        df_gene[cols_finais_gene + outras_gene],
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            col_score: st.column_config.NumberColumn("AlphaMissense Score", format="%.3f"),
                            col_class_dbsnp: st.column_config.TextColumn("Significância Clínica"),
                            "Posicao_Gene": st.column_config.TextColumn("Posição no Gene / cDNA"),
                            "Posicao_Genomica": st.column_config.TextColumn("Posição Genômica / Chr"),
                            "Consequencia_Gene": st.column_config.TextColumn("Consequência Funcional"),
                            "Mudanca_Proteina": st.column_config.TextColumn("Mudança na Proteína")
                        }
                    )

                    st.markdown("---")
                    st.markdown("##### 🧬 Sequência FASTA com Sinalização da Mutação")

                    if fasta_original:
                        linhas_fasta = fasta_original.strip().split("\n")
                        seq_limpa = "".join(linhas_fasta[1:]).replace("\n", "").replace(" ", "")
                        df_mutacoes = df_gene.dropna(subset=['Posicao_Proteina', 'AA_Ref', 'AA_Alt']).copy()

                        if not df_mutacoes.empty:
                            opcoes_mut = [
                                f"{row[col_rs]} | {row['Mudanca_Proteina']} (Posição {int(row['Posicao_Proteina'])})" 
                                for _, row in df_mutacoes.iterrows()
                            ]
                            
                            mut_selecionada = st.selectbox("🎯 Escolha a mutação para destacar no FASTA:", opcoes_mut)
                            idx_sel = opcoes_mut.index(mut_selecionada)
                            row_sel = df_mutacoes.iloc[idx_sel]

                            pos = int(row_sel['Posicao_Proteina']) - 1
                            aa_alt = str(row_sel['AA_Alt'])

                            if 0 <= pos < len(seq_limpa):
                                aa_real_fasta = seq_limpa[pos]
                                st.info(f"📍 **Alteração Selecionada:** Posição `{pos + 1}` | Referência (`{aa_real_fasta}`) ➡️ Alterado (`{aa_alt}`)")

                                seq_com_destaque = (
                                    seq_limpa[:pos] +
                                    f'<span style="background-color: #EF4444; color: #FFFFFF; font-weight: 800; padding: 2px 6px; border-radius: 4px; box-shadow: 0 0 10px rgba(239, 68, 68, 0.8);" title="Mutação: {aa_real_fasta} -> {aa_alt}">{aa_real_fasta} ({aa_alt})</span>' +
                                    seq_limpa[pos + 1:]
                                )

                                blocos = [seq_com_destaque[i:i+60] for i in range(0, len(seq_com_destaque), 60)]
                                fasta_formatado = "<br>".join(blocos)

                                st.markdown(f"""
                                <div style="
                                    background-color: #0F172A; 
                                    border: 1px solid #1E293B; 
                                    border-radius: 12px; 
                                    padding: 16px; 
                                    font-family: 'Fira Code', 'Courier New', monospace; 
                                    font-size: 0.9rem; 
                                    color: #94A3B8; 
                                    line-height: 1.8;
                                    max-height: 350px;
                                    overflow-y: auto;
                                    word-break: break-all;
                                ">
                                    {fasta_formatado}
                                </div>
                                """, unsafe_allow_html=True)
                            else:
                                st.warning("A posição da mutação excede a extensão da sequência FASTA carregada.")
                        else:
                            st.warning("Não há variantes com posições de aminoácidos válidas para destacar.")
                    else:
                        st.warning("Sequência FASTA não disponível no UniProt para esta proteína.")

                with tab_resumo:
                    st.markdown("##### 📋 Resumo dos Atributos e Anotações UniProt")
                    if anotacoes_uniprot:
                        c_res1, c_res2 = st.columns(2)
                        with c_res1:
                            st.markdown("**Domínios & Regiões:**")
                            if anotacoes_uniprot.get("domains"):
                                for d in anotacoes_uniprot["domains"]:
                                    st.markdown(f"- {d}")
                            else:
                                st.write("Nenhum domínio específico listado.")
                                
                            st.markdown("**Sítios de Ligação / Ativos:**")
                            if anotacoes_uniprot.get("binding_sites"):
                                for b in anotacoes_uniprot["binding_sites"]:
                                    st.markdown(f"- {b}")
                            else:
                                st.write("Nenhum sítio ativo listado.")

                        with c_res2:
                            st.markdown("**Peptídeo Sinal:**")
                            if anotacoes_uniprot.get("signal_peptide"):
                                for s in anotacoes_uniprot["signal_peptide"]:
                                    st.markdown(f"- {s}")
                            else:
                                st.write("Sem indicação de peptídeo sinal.")

                            st.markdown("**Interações Proteicas (Rede):**")
                            if anotacoes_uniprot.get("interactions"):
                                for inter in anotacoes_uniprot["interactions"]:
                                    st.markdown(f"- {inter}")
                            else:
                                st.write("Sem interações documentadas.")
                    else:
                        st.info("Anotações detalhadas não disponíveis para este UniProt ID.")

                with tab_3d:
                    st.markdown("##### 🧊 Análise da Estrutura 3D e Rede de Contatos Espaciais (AlphaFold)")
                    pdb_data = obter_pdb_alphafold(uniprot_val)

                    if pdb_data:
                        df_ca = extrair_coordenadas_ca(pdb_data)
                        if not df_ca.empty:
                            G, n_edges, df_lig = calcular_rede_ligacoes(df_ca)
                            
                            m1, m2 = st.columns(2)
                            m1.metric("Resíduos Estruturados (C-Alpha)", len(df_ca))
                            m2.metric("Contatos Espaciais Inter-resíduos (< 8.0Å)", n_edges)
                            
                            st.markdown("**Top Resíduos com Maior Conectividade Espacial (Hubs):**")
                            st.dataframe(df_lig.head(10), use_container_width=True, hide_index=True)
                        else:
                            st.warning("Não foi possível extrair coordenadas dos átomos C-Alpha do PDB.")
                    else:
                        st.warning("Estrutura PDB 3D não encontrada na base do AlphaFold.")