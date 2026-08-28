import os
import pandas as pd

def gerar_interface_html_proteinas():
    tsv_path = os.path.join("resultados", "tabela_comparativa_detalhada_proteinas.tsv")
    html_output_path = os.path.join("resultados", "dashboard_proteinas.html")

    if not os.path.exists(tsv_path):
        print(f"[ERRO] Arquivo {tsv_path} não foi encontrado. Execute o pipeline de dados primeiro.")
        return

    print("Carregando dados do TSV...")
    df = pd.read_csv(tsv_path, sep="\t").fillna("-")

    # Cálculos para os Cards do Dashboard
    total_mutaçoes = len(df)
    total_genes = df['Gene'].nunique() if 'Gene' in df.columns else 0
    total_pdbs = df['PDB_IDs'].apply(lambda x: 0 if x in ['-', 'Sem Estrutura PDB'] else len(x.split(','))).sum()
    
    # Preparar linhas da tabela para HTML
    table_rows_html = []
    for idx, row in df.iterrows():
        gene = row.get('Gene', '-')
        uniprot = row.get('UniProt_ID', '-')
        nome_prot = row.get('Nome_Oficial_Proteina', '-')
        pdb_ids = row.get('PDB_IDs', '-')
        rsid = row.get('rsID_dbSNP', '-')
        troca = row.get('Troca_Aminoacido', '-')
        tipo_mod = row.get('Tipo_Modificacao_FisicoQuimica', '-')
        am_class = row.get('am_classificacao', row.get('am_class', '-'))
        am_score = row.get('am_pathogenicity', row.get('am_score', '-'))
        
        link_uniprot = row.get('Link_UniProt', '-')
        link_rcsb = row.get('Link_RCSB_PDB', '-')
        link_dbsnp = row.get('Link_dbSNP', '-')

        # PDB primário para a visualização 3D
        primeiro_pdb = pdb_ids.split(',')[0].strip() if pdb_ids not in ['-', 'Sem Estrutura PDB'] else ''

        # Badges visuais para o tipo de alteração
        badge_class = "badge-neutral"
        if "Não-Conservativa" in str(tipo_mod):
            badge_class = "badge-danger"
        elif "Conservativa" in str(tipo_mod):
            badge_class = "badge-success"

        # Badges para AlphaMissense
        am_badge = "badge-neutral"
        if "pathogenic" in str(am_class).lower():
            am_badge = "badge-danger"
        elif "benign" in str(am_class).lower():
            am_badge = "badge-success"

        row_html = f"""
        <tr data-pdb="{primeiro_pdb}">
            <td><strong>{gene}</strong></td>
            <td>
                {uniprot}
                {f'<a href="{link_uniprot}" target="_blank" class="external-link" title="Ver no UniProt">↗</a>' if link_uniprot != '-' else ''}
            </td>
            <td><span class="text-truncate" style="max-width: 180px;" title="{nome_prot}">{nome_prot}</span></td>
            <td>
                {rsid}
                {f'<a href="{link_dbsnp}" target="_blank" class="external-link" title="Ver no dbSNP">↗</a>' if link_dbsnp != '-' else ''}
            </td>
            <td><code>{troca}</code></td>
            <td><span class="badge {badge_class}">{tipo_mod}</span></td>
            <td><span class="badge {am_badge}">{am_class}</span> <small>({am_score})</small></td>
            <td>
                {pdb_ids if pdb_ids != '-' else '<span class="text-muted">Ausente</span>'}
                {f'<br><button class="btn-3d" onclick="carregarEstrutura3D(\'{primeiro_pdb}\', \'{gene}\')">👁️ Ver 3D</button>' if primeiro_pdb else ''}
            </td>
        </tr>
        """
        table_rows_html.append(row_html)

    html_content = f"""<!DOCTYPE html>
<html lang="pt-BR">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Painel de Variantes Proteicas & Estruturas 3D</title>
    
    <!-- DataTables CSS & Bootstrap 5 -->
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/1.13.6/css/dataTables.bootstrap5.min.css">
    <link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.4.1/css/buttons.bootstrap5.min.css">

    <!-- 3Dmol.js para Renderização 3D das Proteínas -->
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>

    <style>
        body {{
            background-color: #f4f6f9;
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #333;
        }}
        .header-title {{
            background: linear-gradient(135deg, #1e3c72 0%, #2a5298 100%);
            color: white;
            padding: 2rem 0;
            margin-bottom: 2rem;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        .card-kpi {{
            border: none;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
            transition: transform 0.2s;
        }}
        .card-kpi:hover {{
            transform: translateY(-3px);
        }}
        .kpi-number {{
            font-size: 2rem;
            font-weight: bold;
            color: #1e3c72;
        }}
        .table-container {{
            background: white;
            padding: 1.5rem;
            border-radius: 10px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.08);
        }}
        #viewport3d {{
            width: 100%;
            height: 400px;
            position: relative;
            background-color: #111;
            border-radius: 8px;
            overflow: hidden;
        }}
        .badge-danger {{ background-color: #e74c3c; color: white; }}
        .badge-success {{ background-color: #2ecc71; color: white; }}
        .badge-neutral {{ background-color: #95a5a6; color: white; }}
        .external-link {{
            text-decoration: none;
            margin-left: 3px;
            font-weight: bold;
        }}
        .btn-3d {{
            border: none;
            background: #2a5298;
            color: white;
            font-size: 0.75rem;
            padding: 2px 8px;
            border-radius: 4px;
            margin-top: 4px;
            cursor: pointer;
        }}
        .btn-3d:hover {{
            background: #1e3c72;
        }}
        .text-truncate {{
            display: inline-block;
            white-space: nowrap;
            overflow: hidden;
            text-overflow: ellipsis;
        }}
    </style>
</head>
<body>

<div class="header-title text-center">
    <h2>Plataforma de Análise de Variantes e Estrutura Proteica 3D</h2>
    <p class="lead mb-0">Integração UniProt, dbSNP, AlphaMissense e RCSB PDB</p>
</div>

<div class="container-fluid px-4">
    <!-- Bloco de KPIs -->
    <div class="row mb-4">
        <div class="col-md-4">
            <div class="card card-kpi p-3 text-center">
                <div class="text-muted">Total de Mutações Analisadas</div>
                <div class="kpi-number">{total_mutaçoes}</div>
            </div>
        </div>
        <div class="col-md-4">
            <div class="card card-kpi p-3 text-center">
                <div class="text-muted">Genes Mapeados</div>
                <div class="kpi-number">{total_genes}</div>
            </div>
        </div>
        <div class="col-md-4">
            <div class="card card-kpi p-3 text-center">
                <div class="text-muted">Estruturas PDB Associadas</div>
                <div class="kpi-number">{total_pdbs}</div>
            </div>
        </div>
    </div>

    <div class="row">
        <!-- Tabela Principal -->
        <div class="col-lg-8 mb-4">
            <div class="table-container">
                <h5 class="mb-3">Tabela de Variantes Proteicas</h5>
                <div class="table-responsive">
                    <table id="tabelaProteinas" class="table table-hover table-striped align-middle" style="width:100%">
                        <thead class="table-dark">
                            <tr>
                                <th>Gene</th>
                                <th>UniProt</th>
                                <th>Proteína</th>
                                <th>rsID</th>
                                <th>Mutação</th>
                                <th>Efeito Físico-Químico</th>
                                <th>AlphaMissense</th>
                                <th>PDB (3D)</th>
                            </tr>
                        </thead>
                        <tbody>
                            {"".join(table_rows_html)}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>

        <!-- Painel de Visualização 3D -->
        <div class="col-lg-4 mb-4">
            <div class="table-container sticky-top" style="top: 20px;">
                <h5 class="mb-2">Visualizador Estrutural 3D</h5>
                <p id="infoPDB" class="text-muted small">Clique no botão "👁️ Ver 3D" de qualquer mutação para carregar a estrutura macromolecular.</p>
                
                <div id="viewport3d"></div>
                
                <div class="mt-3 text-end">
                    <button class="btn btn-sm btn-outline-secondary" onclick="resetarVisao3D()">Resetar Vista</button>
                </div>
            </div>
        </div>
    </div>
</div>

<!-- Scripts JavaScript para Interatividade -->
<script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/js/bootstrap.bundle.min.js"></script>

<!-- DataTables + Plugins de Exportação -->
<script src="https://cdn.datatables.net/1.13.6/js/jquery.dataTables.min.js"></script>
<script src="https://cdn.datatables.net/1.13.6/js/dataTables.bootstrap5.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.1/js/dataTables.buttons.min.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.1/js/buttons.bootstrap5.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/pdfmake.min.js"></script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/pdfmake/0.1.53/vfs_fonts.js"></script>
<script src="https://cdn.datatables.net/buttons/2.4.1/js/buttons.html5.min.js"></script>

<script>
    let viewer3D = null;

    $(document).ready(function() {{
        // Inicializa DataTables com busca, ordenação e exportação
        $('#tabelaProteinas').DataTable({{
            language: {{
                url: 'https://cdn.datatables.net/plug-ins/1.13.6/i18n/pt-BR.json'
            }},
            dom: 'Bfrtip',
            buttons: [
                {{ extend: 'csv', className: 'btn btn-sm btn-outline-primary me-1' }},
                {{ extend: 'excel', className: 'btn btn-sm btn-outline-success me-1' }},
                {{ extend: 'pdf', className: 'btn btn-sm btn-outline-danger me-1' }}
            ],
            pageLength: 10,
            order: [[0, 'asc']]
        }});

        // Inicializa o visualizador 3Dmol
        let element = $('#viewport3d');
        let config = {{ backgroundColor: 'black' }};
        viewer3D = $3Dmol.createViewer(element, config);
    }});

    // Função para buscar a estrutura PDB diretamente do servidor do RCSB PDB e renderizar em 3D
    function carregarEstrutura3D(pdbId, gene) {{
        if (!pdbId) return;

        $('#infoPDB').html('Carregando estrutura PDB: <strong>' + pdbId + '</strong> (' + gene + ')...');

        let pdbUrl = 'https://files.rcsb.org/download/' + pdbId + '.pdb';

        viewer3D.clear();
        $.get(pdbUrl, function(data) {{
            viewer3D.addModel(data, "pdb");
            viewer3D.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
            viewer3D.zoomTo();
            viewer3D.render();
            $('#infoPDB').html('Exibindo PDB: <strong>' + pdbId + '</strong> | Gene: <strong>' + gene + '</strong>');
        }}).fail(function() {{
            $('#infoPDB').html('<span class="text-danger">Erro ao carregar o arquivo PDB ' + pdbId + '.</span>');
        }});
    }}

    function resetarVisao3D() {{
        if (viewer3D) {{
            viewer3D.zoomTo();
            viewer3D.render();
        }}
    }}
</script>

</body>
</html>
"""

    with open(html_output_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"\n[SUCESSO] Interface HTML criada com sucesso em: '{html_output_path}'!")

if __name__ == "__main__":
    gerar_interface_html_proteinas()