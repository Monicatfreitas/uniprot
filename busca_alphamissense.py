cat << 'EOF' > gerar_interface.py
import streamlit as st
import pandas as pd
import numpy as np

st.set_page_config(page_title="Consolidado por Gene", layout="wide")

@st.cache_data
def carregar_dados():
    caminhos = [
        "resultados/tabela_correlacao_genes_uniprot.tsv",
        "resultados/tabela_comparativa_detalhada_proteinas.tsv",
        "resultados/dbsnp_alphamissense_unificado.tsv"
    ]
    df = pd.DataFrame()
    for caminho in caminhos:
        try:
            df = pd.read_csv(caminho, sep="\t", on_bad_lines="skip")
            if not df.empty:
                break
        except Exception:
            continue
    df.columns = df.columns.str.strip()
    return df

df = carregar_dados()

st.title("🧬 Resumo Compilado por Gene (dbSNP + AlphaMissense)")

if df.empty:
    st.error("Nenhum arquivo de dados foi encontrado na pasta 'resultados/'.")
else:
    col_gene = next((c for c in df.columns if c.lower() in ['gene', 'gene_name', 'gene name', 'symbol']), None)

    if not col_gene:
        st.error("Coluna de Gene não encontrada no arquivo.")
    else:
        gene_busca = st.text_input("Digite o nome do Gene:", placeholder="Ex: IL17A, RIT1, TP53...").strip()

        if gene_busca:
            df_gene = df[df[col_gene].astype(str).str.upper() == gene_busca.upper()].copy()

            if not df_gene.empty:
                st.success(f"Gene **{gene_busca.upper()}** localizado com sucesso!")

                def limpar_e_converter(serie):
                    s = serie.astype(str).str.replace(',', '.').str.strip()
                    s = s.replace(['nan', 'none', 'n/a', '-', 'null', '', '<na>'], np.nan)
                    return pd.to_numeric(s, errors='coerce')

                # CARDS DE MÉTRICAS
                c1, c2, c3, c4 = st.columns(4)
                with c1:
                    st.metric("Total de Mutações", len(df_gene))

                with c2:
                    col_rs = next((c for c in df.columns if 'rs' in c.lower() or 'dbsnp' in c.lower()), None)
                    if col_rs:
                        rs_count = df_gene[col_rs].dropna().astype(str).str.strip()
                        rs_count = rs_count[~rs_count.isin(['-', 'nan', '', 'N/A'])].nunique()
                        st.metric("Variantes dbSNP (rsID)", rs_count)
                    else:
                        st.metric("Variantes dbSNP", "N/A")

                with c3:
                    col_freq = next((c for c in df.columns if any(k in c.lower() for k in ['freq', 'af', 'allele_freq', 'frequencia'])), None)
                    if col_freq:
                        vals_f = limpar_e_converter(df_gene[col_freq])
                        m_freq = vals_f.mean()
                        st.metric("Frequência Média", f"{m_freq:.4f}" if not np.isnan(m_freq) else "N/A")
                    else:
                        st.metric("Frequência Média", "N/A")

                with c4:
                    col_score = next((c for c in df.columns if any(k in c.lower() for k in ['am_pathogenicity', 'alphamissense', 'score'])), None)
                    if col_score:
                        vals_s = limpar_e_converter(df_gene[col_score])
                        m_score = vals_s.mean()
                        st.metric("Score Médio AlphaMissense", f"{m_score:.3f}" if not np.isnan(m_score) else "N/A")
                    else:
                        st.metric("Score AlphaMissense", "N/A")

                st.divider()

                # RESUMO COMPILADO
                st.subheader("📋 Resumo dos Campos e Atributos")
                resumo_list = []
                for col in df_gene.columns:
                    serie = df_gene[col].astype(str).str.strip()
                    validos = [str(x) for x in serie[~serie.isin(['-', 'nan', 'None', '', 'N/A', '<NA>'])].unique()]
                    num_vals = limpar_e_converter(df_gene[col]).dropna()

                    if len(num_vals) > 0 and len(num_vals) == len(validos):
                        info = f"Média: {num_vals.mean():.4f} | Mín: {num_vals.min()} | Máx: {num_vals.max()}"
                    elif len(validos) > 0:
                        if len(validos) <= 5:
                            info = ", ".join(validos)
                        else:
                            info = f"{len(validos)} valores únicos (Ex: {', '.join(validos[:5])}...)"
                    else:
                        info = "Sem dados registrados para este gene"

                    resumo_list.append({"Coluna / Informação": col, "Resumo Compilado": info})

                st.table(pd.DataFrame(resumo_list))

                st.divider()

                # TABELA DETALHADA
                st.subheader("🔍 Tabela Detalhada com Todas as Mutações Registradas")
                st.dataframe(df_gene, use_container_width=True)

            else:
                st.warning(f"Nenhum registro encontrado para o gene '{gene_busca}'.")
        else:
            st.info("Digite o nome de um Gene acima para carregar as informações.")
EOF