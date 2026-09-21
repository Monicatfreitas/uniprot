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

# --- ESTILIZAÇÃO CUSTOMIZADA (CSS) ---
st.markdown("""
<style>
    /* Reset e Tipografia do Sistema */
    html, body, [class*="css"] {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    }
    
    /* Título Principal da App */
    .app-title {
        font-size: 1.8rem;
        font-weight: 700;
        color: #F8FAFC;
        letter-spacing: -0.02em;
        margin-bottom: 1.5rem;
        padding-bottom: 0.5rem;
        border-bottom: 1px solid #1E293B;
    }
    
    /* Card de Cabeçalho do Gene */
    .header-card {
        background-color: #0F172A;
        border: 1px solid #1E293B;
        border-radius: 8px;
        padding: 20px 24px;
        margin-bottom: 24px;
    }
    .gene-label {
        font-size: 0.85rem;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        color: #64748B;
        font-weight: 600;
    }
    .gene-title-container {
        display: flex;
        align-items: center;
        gap: 12px;
        margin-top: 4px;
    }
    .gene-title {
        font-size: 1.75rem;
        font-weight: 700;
        color: #F8FAFC;
        margin: 0;
    }
    .uniprot-badge {
        background-color: #1E293B;
        color: #38BDF8;
        border: 1px solid #334155;
        padding: 4px 12px;
        border-radius: 4px;
        font-family: "JetBrains Mono", monospace;
        font-size: 0.85rem;
        font-weight: 600;
    }
    .protein-subtitle {
        color: #94A3B8;
        font-size: 0.95rem;
        margin-top: 8px;
    }

    /* Cards de Métricas (st.metric) */
    div[data-testid="stMetric"] {
        background-color: #0F172A;
        border: 1px solid #1E293B;
        padding: 16px;
        border-radius: 8px;
    }
    div[data-testid="stMetricLabel"] {
        color: #64748B !important;
        font-size: 0.8rem !important;
        font-weight: 600 !important;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    div[data-testid="stMetricValue"] {
        color: #F8FAFC !important;
        font-size: 1.6rem !important;
        font-weight: 700 !important;
    }
    
    /* Abas de Navegação (st.tabs) */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
        border-bottom: 1px solid #1E293B;
    }
    .stTabs [data-baseweb="tab"] {
        padding: 8px 16px;
        border-radius: 4px 4px 0 0;
        font-weight: 500;
        font-size: 0.9rem;
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

st.markdown('<div class="app-title">Painel de Análise Genômica e Estrutural</div>', unsafe_allow_html=True)

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
                    <div class="gene-label">Gene Selecionado</div>
                    <div class="gene-title-container">
                        <h1 class="gene-title">{gene_selecionado}</h1>
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
                st.subheader(" Painel Visual da Proteína")
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
                    hotspots = contagem_sitio[contagem_sitio > 1].sort_values(ascending=False)

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
                            zerolinecolor='#4A5568',
                            
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
                st.subheader(" Detalhamento das Alterações")
                tab_proteina, tab_gene, tab_resumo, tab_3d = st.tabs([
                    " Alterações na Proteína",
                    " Gene, Genômica & FASTA",
                    " Atributos & Anotações",
                    " Estrutura 3D, Sinalização de Peptídeo Sinal e Rede de Contatos Espaciais"
                ])

                colunas_redundantes = ['AlphaMissense_Classificacao', 'AM_Score_Num']

                with tab_proteina:
                    st.markdown("#####  Mutações com Foco na Estrutura Proteica e Aminoácidos")
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
                    st.markdown("#####  Mutações com Foco em Posição no Gene (cDNA/CDS), Genômica e Consequência")
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
                    st.markdown("#####  Sequência FASTA com Sinalização da Mutação")

                    if fasta_original:
                        linhas_fasta = fasta_original.strip().split("\n")
                        seq_limpa = "".join(linhas_fasta[1:]).replace("\n", "").replace(" ", "")
                        df_mutacoes = df_gene.dropna(subset=['Posicao_Proteina', 'AA_Ref', 'AA_Alt']).copy()

                        if not df_mutacoes.empty:
                            opcoes_mut = [
                                f"{row[col_rs]} | {row['Mudanca_Proteina']} (Posição {int(row['Posicao_Proteina'])})" 
                                for _, row in df_mutacoes.iterrows()
                            ]
                            
                            mut_selecionada = st.selectbox(" Escolha a mutação para destacar no FASTA:", opcoes_mut)
                            idx_sel = opcoes_mut.index(mut_selecionada)
                            row_sel = df_mutacoes.iloc[idx_sel]

                            pos = int(row_sel['Posicao_Proteina']) - 1
                            aa_alt = str(row_sel['AA_Alt'])

                            if 0 <= pos < len(seq_limpa):
                                aa_real_fasta = seq_limpa[pos]
                                st.info(f"📍 **Alteração Selecionada:** Posição `{pos + 1}` | Referência (`{aa_real_fasta}`) ➡️ Alterado (`{aa_alt}`)")

                with tab_resumo:
                    st.markdown("#####  Resumo dos Atributos & Anotações Funcionais (UniProt)")
                    
                    if anotacoes_uniprot:
                        col_a1, col_a2 = st.columns(2)
                        with col_a1:
                            st.markdown("** Peptídeo Sinal:**")
                            if anotacoes_uniprot.get("signal_peptide"):
                                for item in anotacoes_uniprot["signal_peptide"]:
                                    st.write(f"- {item}")
                            else:
                                st.write("Nenhum peptídeo sinal mapeado.")

                            st.markdown("** Domínios e Regiões:**")
                            if anotacoes_uniprot.get("domains"):
                                for item in anotacoes_uniprot["domains"]:
                                    st.write(f"- {item}")
                            else:
                                st.write("Nenhum domínio específico encontrado.")

                        with col_a2:
                            st.markdown("** Sítios de Ligação / Ativos:**")
                            if anotacoes_uniprot.get("binding_sites"):
                                for item in anotacoes_uniprot["binding_sites"]:
                                    st.write(f"- {item}")
                            else:
                                st.write("Nenhum sítio de ligação informado.")

                            st.markdown("**🔗 Interações Proteína-Proteína:**")
                            if anotacoes_uniprot.get("interactions"):
                                for item in set(anotacoes_uniprot["interactions"]):
                                    st.write(f"- {item}")
                            else:
                                st.write("Nenhuma interação registrada.")
                    else:
                        st.info("Anotações da UniProt não disponíveis para este gene.")

                with tab_3d:
                    if uniprot_val == "N/A":
                        st.info("UniProt ID não informado para carregar os dados de estrutura.")
                    else:
                        pdb_code = obter_pdb_alphafold(uniprot_val)

                        if pdb_code:
                            df_ca = extrair_coordenadas_ca(pdb_code)
                            G, total_edges, df_lig = calcular_rede_ligacoes(df_ca, distancia_corte=8.0)

                            col3d_1, col3d_2 = st.columns([1, 1])

                            with col3d_1:
                                st.markdown("######  Estrutura 3D / Backbone C-Alpha (AlphaFold)")
                                
                                sig_range = anotacoes_uniprot.get("signal_range")
                                if sig_range:
                                    st.caption(f"🔴 **Peptídeo Sinal nos Resíduos {sig_range[0]}-{sig_range[1]}**")
                                else:
                                    st.caption("ℹ️ Nenhum peptídeo sinal mapeado para esta sequência.")

                                if not df_ca.empty:
                                    # Plotly 3D Backbone
                                    fig_3d = go.Figure()
                                    
                                    # Conexão contínua
                                    fig_3d.add_trace(go.Scatter3d(
                                        x=df_ca['x'], y=df_ca['y'], z=df_ca['z'],
                                        mode='lines',
                                        line=dict(color='#64748B', width=2),
                                        hoverinfo='none',
                                        showlegend=False
                                    ))

                                    # Se houver peptídeo sinal
                                    if sig_range:
                                        df_sig = df_ca[(df_ca['res_num'] >= sig_range[0]) & (df_ca['res_num'] <= sig_range[1])]
                                        df_rest = df_ca[(df_ca['res_num'] < sig_range[0]) | (df_ca['res_num'] > sig_range[1])]
                                        
                                        if not df_sig.empty:
                                            fig_3d.add_trace(go.Scatter3d(
                                                x=df_sig['x'], y=df_sig['y'], z=df_sig['z'],
                                                mode='markers',
                                                marker=dict(size=4, color='#EF4444'),
                                                name='Peptídeo Sinal',
                                                text=[f"{r['res_name']}{r['res_num']}" for _, r in df_sig.iterrows()],
                                                hovertemplate="<b>Peptídeo Sinal</b><br>Resíduo: %{text}<extra></extra>"
                                            ))
                                        if not df_rest.empty:
                                            fig_3d.add_trace(go.Scatter3d(
                                                x=df_rest['x'], y=df_rest['y'], z=df_rest['z'],
                                                mode='markers',
                                                marker=dict(size=3, color='#3B82F6'),
                                                name='Cadeia Proteica',
                                                text=[f"{r['res_name']}{r['res_num']}" for _, r in df_rest.iterrows()],
                                                hovertemplate="<b>Resíduo:</b> %{text}<extra></extra>"
                                            ))
                                    else:
                                        fig_3d.add_trace(go.Scatter3d(
                                            x=df_ca['x'], y=df_ca['y'], z=df_ca['z'],
                                            mode='markers',
                                            marker=dict(size=3, color='#3B82F6'),
                                            name='Resíduos C-Alpha',
                                            text=[f"{r['res_name']}{r['res_num']}" for _, r in df_ca.iterrows()],
                                            hovertemplate="<b>Resíduo:</b> %{text}<extra></extra>"
                                        ))

                                    fig_3d.update_layout(
                                        height=480,
                                        margin=dict(l=0, r=0, t=0, b=0),
                                        scene=dict(
                                            xaxis=dict(visible=False),
                                            yaxis=dict(visible=False),
                                            zaxis=dict(visible=False),
                                            bgcolor='rgba(0,0,0,0)'
                                        )
                                    )
                                    st.plotly_chart(fig_3d, use_container_width=True)

                            with col3d_2:
                                st.markdown("###### 🌐 Rede de Ligações Espaciais (C-Alpha < 8Å)")

                                mc1, mc2 = st.columns(2)
                                with mc1:
                                    st.metric("Total de Resíduos Mapeados", f"{len(df_ca):,}".replace(",", "."))
                                with mc2:
                                    st.metric("Total de Contatos/Ligações Spatial", f"{total_edges:,}".replace(",", "."))

                                if not df_ca.empty and G.number_of_edges() > 0:
                                    # Rede Tridimensional 3D
                                    edge_x, edge_y, edge_z = [], [], []
                                    pos_dict = {row['res_num']: (row['x'], row['y'], row['z']) for _, row in df_ca.iterrows()}

                                    for edge in G.edges():
                                        x0, y0, z0 = pos_dict[edge[0]]
                                        x1, y1, z1 = pos_dict[edge[1]]
                                        edge_x.extend([x0, x1, None])
                                        edge_y.extend([y0, y1, None])
                                        edge_z.extend([z0, z1, None])

                                    edge_trace = go.Scatter3d(
                                        x=edge_x, y=edge_y, z=edge_z,
                                        line=dict(width=1, color='rgba(150, 150, 150, 0.3)'),
                                        hoverinfo='none', mode='lines', showlegend=False
                                    )

                                    node_x, node_y, node_z, node_degrees, node_text = [], [], [], [], []
                                    for node in G.nodes():
                                        x, y, z = pos_dict[node]
                                        node_x.append(x)
                                        node_y.append(y)
                                        node_z.append(z)
                                        deg = G.degree(node)
                                        node_degrees.append(deg)
                                        node_text.append(f"{G.nodes[node]['name']} ({deg} ligações 3D)")

                                    node_trace = go.Scatter3d(
                                        x=node_x, y=node_y, z=node_z,
                                        mode='markers',
                                        hoverinfo='text',
                                        text=node_text,
                                        marker=dict(
                                            showscale=True,
                                            colorscale='Plasma',
                                            color=node_degrees,
                                            size=5,
                                            colorbar=dict(thickness=12, title='Ligações 3D', x=1.0)
                                        ),
                                        showlegend=False
                                    )

                                    fig_net = go.Figure(data=[edge_trace, node_trace])
                                    fig_net.update_layout(
                                        height=480,
                                        margin=dict(l=0, r=0, t=0, b=0),
                                        scene=dict(
                                            xaxis=dict(visible=False),
                                            yaxis=dict(visible=False),
                                            zaxis=dict(visible=False),
                                            bgcolor='rgba(0,0,0,0)'
                                        )
                                    )
                                    st.plotly_chart(fig_net, use_container_width=True)

                                st.markdown("###### Resíduos com Maior Conectividade (*Hubs*):")
                                if not df_lig.empty:
                                    st.dataframe(
                                        df_lig.head(10),
                                        use_container_width=True,
                                        hide_index=True
                                    )
                        else:
                            st.error("Não foi possível carregar a estrutura PDB correspondente no AlphaFold DB.")