import json
import os
import pandas as pd


def gerar_exportacao_tabelas_interface(
    tsv_input="resultados/tabela_comparativa_detalhada_proteinas.tsv",
    json_output="resultados/dados_genes_tabela.json",
    html_output="resultados/visualizacao_interface.html",
):
    if not os.path.exists(tsv_input):
        print(f"[ERRO] Arquivo de entrada não encontrado em: {tsv_input}")
        return

    print("Lendo tabela comparativa unificada...")
    df = pd.read_csv(tsv_input, sep="\t")
    df = df.fillna("-")

    estrutura_tabelas_genes = []

    print("Agrupando mutações e scores por gene...")
    for gene, group in df.groupby("Gene"):
        primeira_linha = group.iloc[0]

        mutacoes_tabela = []
        for _, row in group.iterrows():
            mutacao = {
                "rsid": str(row.get("rsID_dbSNP", "-")),
                "proteina_hgvs": str(row.get("Proteina_HGVS", "-")),
                "troca_aminoacido": str(row.get("Troca_Aminoacido", "-")),
                "aa_referencia": str(row.get("AA_Referencia", "-")),
                "aa_alterado": str(row.get("AA_Alterado", "-")),
                "posicao": str(row.get("Posicao_Proteina", "-")),
                "tipo_modificacao": str(
                    row.get("Tipo_Modificacao_FisicoQuimica", "-")
                ),
                "comparativo_fontes": str(row.get("Comparativo_Fontes", "-")),
                "link_dbsnp": str(row.get("Link_dbSNP", "-")),
                "alphamissense_score": str(
                    row.get(
                        "am_pathogenicity", row.get("pathogenicity_score", "-")
                    )
                ),
                "alphamissense_class": str(
                    row.get("am_class", row.get("class", "-"))
                ),
            }
            mutacoes_tabela.append(mutacao)

        gene_obj = {
            "gene": str(gene),
            "uniprot_id": str(primeira_linha.get("UniProt_ID", "-")),
            "nome_proteina": str(
                primeira_linha.get("Nome_Oficial_Proteina", "-")
            ),
            "link_uniprot": str(primeira_linha.get("Link_UniProt", "-")),
            "total_mutacoes": len(mutacoes_tabela),
            "tabela_mutacoes": mutacoes_tabela,
        }
        estrutura_tabelas_genes.append(gene_obj)

    # Exporta JSON
    with open(json_output, "w", encoding="utf-8") as f:
        json.dump(estrutura_tabelas_genes, f, indent=2, ensure_ascii=False)

    # Exporta HTML
    html_content = """<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <title>Painel Comparativo de Mutações por Gene</title>
    <style>
        body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; margin: 30px; background-color: #f4f6f8; color: #333; }
        h1 { color: #1a252f; margin-bottom: 25px; }
        .gene-card { background: white; border-radius: 8px; padding: 20px; margin-bottom: 30px; box-shadow: 0 4px 6px rgba(0,0,0,0.05); border: 1px solid #e1e8ed; }
        .gene-header { border-bottom: 2px solid #3498db; padding-bottom: 12px; margin-bottom: 15px; }
        .gene-header h2 { margin: 0 0 8px 0; color: #2c3e50; }
        .gene-header p { margin: 0; color: #576574; font-size: 14px; }
        table { width: 100%; border-collapse: collapse; margin-top: 12px; font-size: 13px; }
        th, td { border: 1px solid #e1e8ed; padding: 8px 12px; text-align: left; }
        th { background-color: #3498db; color: white; font-weight: 600; }
        tr:nth-child(even) { background-color: #f8f9fa; }
        tr:hover { background-color: #f1f2f6; }
        a { color: #2980b9; text-decoration: none; font-weight: bold; }
        a:hover { text-decoration: underline; }
        .badge { display: inline-block; padding: 3px 7px; border-radius: 4px; font-size: 11px; font-weight: 600; background: #e0e0e0; color: #333; }
    </style>
</head>
<body>
    <h1>Painel de Genes e Mutações Proteicas</h1>
"""

    for g in estrutura_tabelas_genes:
        html_content += f"""
    <div class="gene-card">
        <div class="gene-header">
            <h2>Gene: {g['gene']}</h2>
            <p><strong>Proteína:</strong> {g['nome_proteina']} | 
               <strong>UniProt ID:</strong> <a href="{g['link_uniprot']}" target="_blank">{g['uniprot_id']}</a> | 
               <strong>Total de Mutações:</strong> {g['total_mutacoes']}</p>
        </div>
        <table>
            <thead>
                <tr>
                    <th>rsID</th>
                    <th>HGVS Proteína</th>
                    <th>Troca AA</th>
                    <th>AA Ref</th>
                    <th>AA Alt</th>
                    <th>Posição</th>
                    <th>Modificação Fisico-Química</th>
                    <th>Fonte</th>
                    <th>Score AM</th>
                    <th>Classe AM</th>
                    <th>Link dbSNP</th>
                </tr>
            </thead>
            <tbody>
"""
        for m in g["tabela_mutacoes"]:
            html_content += f"""
                <tr>
                    <td><strong>{m['rsid']}</strong></td>
                    <td>{m['proteina_hgvs']}</td>
                    <td>{m['troca_aminoacido']}</td>
                    <td>{m['aa_referencia']}</td>
                    <td>{m['aa_alterado']}</td>
                    <td>{m['posicao']}</td>
                    <td>{m['tipo_modificacao']}</td>
                    <td><span class="badge">{m['comparativo_fontes']}</span></td>
                    <td>{m['alphamissense_score']}</td>
                    <td>{m['alphamissense_class']}</td>
                    <td><a href="{m['link_dbsnp']}" target="_blank">NCBI</a></td>
                </tr>
"""
        html_content += """
            </tbody>
        </table>
    </div>
"""

    html_content += "</body></html>"

    with open(html_output, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"\n[SUCESSO] Arquivo JSON gerado em: '{json_output}'")
    print(f"[SUCESSO] Relatório HTML gerado em: '{html_output}'")


if __name__ == "__main__":
    gerar_exportacao_tabelas_interface()