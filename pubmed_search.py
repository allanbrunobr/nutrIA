import os
import streamlit as st
from Bio import Entrez
from typing import List, Dict
import pandas as pd
from datetime import datetime
import time
from openai import OpenAI
import xml.etree.ElementTree as ET
import logging

# Configure logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler('pubmed_app.log')
    ]
)
logger = logging.getLogger(__name__)


class PubMedSearch:
    def __init__(self, email: str, pubmed_key: str, openai_key: str = None):
        self.email = email
        self.pubmed_key = pubmed_key
        self.openai_key = openai_key
        Entrez.email = email
        Entrez.api_key = pubmed_key
        if openai_key:
            self.openai_client = OpenAI(api_key=openai_key)
            logger.debug("OpenAI client initialized successfully")
        else:
            logger.warning("OpenAI key not provided")

    def search_articles(self,
                        query: str,
                        max_results: int = 10,
                        start_year: int = None,
                        end_year: int = None) -> List[Dict]:
        """Realiza busca no PubMed"""
        logger.debug(f"Starting search with query: {query}")

        if start_year and end_year:
            query = f"{query} AND ({start_year}:{end_year}[pdat])"
            logger.debug(f"Modified query with year range: {query}")

        try:
            with st.spinner('Buscando artigos...'):
                logger.debug("Executing PubMed search")
                search_handle = Entrez.esearch(
                    db="pubmed",
                    term=query,
                    retmax=max_results,
                    sort="relevance"
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()

                if not search_results["IdList"]:
                    logger.debug("No results found")
                    return []

                logger.debug(f"Found {len(search_results['IdList'])} articles")
                fetch_handle = Entrez.efetch(
                    db="pubmed",
                    id=search_results["IdList"],
                    rettype="abstract",
                    retmode="xml"
                )

                articles = Entrez.read(fetch_handle)['PubmedArticle']
                fetch_handle.close()

                processed_articles = []
                for article in articles:
                    try:
                        article_data = article['MedlineCitation']['Article']

                        # Extrai autores
                        authors = []
                        if 'AuthorList' in article_data:
                            for author in article_data['AuthorList']:
                                if 'LastName' in author and 'ForeName' in author:
                                    authors.append(f"{author['LastName']} {author['ForeName']}")

                        # Extrai abstract
                        abstract = ""
                        if 'Abstract' in article_data:
                            abstract_text = article_data['Abstract']['AbstractText']
                            if isinstance(abstract_text, list):
                                abstract = ' '.join(str(text) for text in abstract_text)
                            else:
                                abstract = str(abstract_text)

                        # Verifica se é artigo gratuito
                        is_free = False
                        pmc_id = None
                        if 'ArticleIdList' in article['PubmedData']:
                            for id_item in article['PubmedData']['ArticleIdList']:
                                if id_item.attributes.get('IdType') == 'pmc':
                                    is_free = True
                                    pmc_id = str(id_item)

                        # Extrai DOI
                        doi = "N/A"
                        if 'ELocationID' in article_data:
                            for id in article_data['ELocationID']:
                                if hasattr(id, 'attributes') and id.attributes.get('EIdType') == 'doi':
                                    doi = str(id)

                        processed_article = {
                            'title': article_data['ArticleTitle'],
                            'authors': '; '.join(authors[:3]) + (' et al.' if len(authors) > 3 else ''),
                            'journal': article_data['Journal']['Title'],
                            'year': int(article_data['Journal']['JournalIssue']['PubDate'].get('Year', 0)),
                            'abstract': abstract[:300] + '...' if len(abstract) > 300 else abstract,
                            'pmid': str(article['MedlineCitation']['PMID']),
                            'doi': doi,
                            'url': f"https://pubmed.ncbi.nlm.nih.gov/{article['MedlineCitation']['PMID']}",
                            'is_free': is_free,
                            'pmc_id': pmc_id
                        }

                        processed_articles.append(processed_article)
                        logger.debug(f"Processed article PMID: {processed_article['pmid']}")

                    except Exception as e:
                        logger.error(f"Error processing article: {str(e)}", exc_info=True)
                        st.error(f"Erro ao processar artigo: {str(e)}")
                        continue

                logger.debug(f"Successfully processed {len(processed_articles)} articles")
                return processed_articles

        except Exception as e:
            logger.error(f"Error in search: {str(e)}", exc_info=True)
            st.error(f"Erro na busca: {str(e)}")
            return []

    def get_full_article(self, pmc_id: str) -> str:
        """Obtém o texto completo do artigo do PMC"""
        logger.debug(f"Iniciando get_full_article para PMC ID: {pmc_id}")
        print(f"Iniciando get_full_article para PMC ID: {pmc_id}")  # Debug print

        if not pmc_id:
            logger.error("PMC ID está vazio ou None")
            print("PMC ID está vazio ou None")  # Debug print
            return ""

        # Remove o prefixo "PMC" se existir
        pmc_id = pmc_id.replace("PMC", "")

        try:
            logger.debug(f"Tentando efetch no PMC com ID: {pmc_id}")
            print(f"Tentando efetch no PMC com ID: {pmc_id}")  # Debug print

            # Primeiro, verifica se o artigo está disponível
            handle = Entrez.efetch(
                db="pmc",
                id=pmc_id,
                rettype="medline",
                retmode="text"
            )

            medline_data = handle.read()
            handle.close()

            if "PMC - Not Found" in medline_data or "PMC - Bad ID" in medline_data:
                logger.error(f"Artigo PMC {pmc_id} não encontrado ou acesso negado")
                print(f"Artigo PMC {pmc_id} não encontrado ou acesso negado")
                st.error("Este artigo não está disponível no PMC ou requer acesso institucional.")
                return ""

            # Se o artigo está disponível, busca o texto completo
            logger.debug("Buscando texto completo em XML")
            print("Buscando texto completo em XML")  # Debug print

            handle = Entrez.efetch(
                db="pmc",
                id=pmc_id,
                rettype="xml",
                retmode="xml"
            )

            xml_content = handle.read()
            handle.close()

            if not xml_content:
                logger.error("XML vazio retornado")
                print("XML vazio retornado")
                st.error("Não foi possível recuperar o conteúdo do artigo.")
                return ""

            root = ET.fromstring(xml_content)

            # Log do XML para debug
            logger.debug(f"XML recebido: {xml_content[:500]}...")
            print(f"XML recebido: {xml_content[:500]}...")

            article_text = []

            # Extract article-meta data
            article_meta = root.find(".//article-meta")
            if article_meta is not None:
                # Get title
                title = article_meta.find(".//article-title")
                if title is not None and title.text:
                    article_text.append(f"TÍTULO: {title.text}")

                # Get abstract
                abstract = article_meta.find(".//abstract")
                if abstract is not None:
                    article_text.append("\nRESUMO:")
                    for p in abstract.findall(".//p"):
                        if p is not None and p.text:
                            article_text.append(p.text)

            # Extract body content
            body = root.find(".//body")
            if body is not None:
                article_text.append("\nCONTEÚDO PRINCIPAL:")

                # Extrai seções
                for sec in body.findall(".//sec"):
                    title = sec.find("title")
                    if title is not None and title.text:
                        article_text.append(f"\n{title.text}:")

                    # Extrai parágrafos da seção
                    for p in sec.findall(".//p"):
                        if p is not None:
                            # Pega todo o texto do parágrafo, incluindo subelementos
                            text = ''.join(p.itertext()).strip()
                            if text:
                                article_text.append(text)

            full_text = "\n".join(article_text)

            if not full_text:
                logger.error("Nenhum texto extraído do artigo")
                print("Nenhum texto extraído do artigo")
                st.error("Não foi possível extrair o texto do artigo.")
                return ""

            logger.info(f"Texto extraído com sucesso. Tamanho: {len(full_text)}")
            print(f"Texto extraído com sucesso. Tamanho: {len(full_text)}")
            return full_text

        except Exception as e:
            error_msg = f"Erro ao extrair conteúdo do PMC: {str(e)}"
            logger.error(error_msg, exc_info=True)
            print(error_msg)
            st.error(error_msg)
            return ""

    def analyze_article(self, text: str, analysis_type: str = "complete") -> str:
        """Analisa o artigo usando OpenAI"""
        logger.debug(f"Starting article analysis with type: {analysis_type}")
        print(f"Starting article analysis with type: {analysis_type}")  # Debug print

        if not self.openai_key:
            error_msg = "Chave da OpenAI não configurada ou inválida"
            logger.error(error_msg)
            print(error_msg)  # Debug print
            st.error(error_msg)
            return None

        if not text:
            error_msg = "Texto do artigo está vazio"
            logger.error(error_msg)
            print(error_msg)  # Debug print
            st.error(error_msg)
            return None

        try:
            # Validar a chave OpenAI
            if not self.openai_key.startswith('sk-') or '--' in self.openai_key:
                error_msg = "Formato da chave OpenAI parece inválido"
                logger.error(error_msg)
                print(error_msg)  # Debug print
                st.error(error_msg)
                return None

            logger.debug("Preparing OpenAI request")
            print("Preparing OpenAI request")  # Debug print

            prompts = {
                "complete": """
                    Analise este artigo científico e forneça:
                    1. Principais descobertas e conclusões
                    2. Metodologia utilizada
                    3. Limitações do estudo
                    4. Implicações práticas
                    5. Recomendações baseadas nas evidências
                """,
                "methodology": "Analise detalhadamente a metodologia usada neste estudo.",
                "results": "Resuma os principais resultados e descobertas do estudo.",
                "practical": """
                    Foque nas implicações práticas:
                    1. Recomendações práticas
                    2. Como aplicar no dia a dia
                    3. Pontos de atenção
                """
            }

            # Limitar o tamanho do texto para evitar exceder limites da API
            max_text_length = 14000  # Aproximadamente 4000 tokens
            text = text[:max_text_length] + "..." if len(text) > max_text_length else text

            logger.debug("Sending request to OpenAI")
            print("Sending request to OpenAI")  # Debug print
            st.info("Enviando requisição para análise...")  # Feedback visual adicional

            try:
                response = self.openai_client.chat.completions.create(
                    model="gpt-4-turbo-preview",
                    messages=[
                        {"role": "system",
                         "content": "Você é um assistente especializado em análise de artigos científicos da área de nutrição."},
                        {"role": "user", "content": f"{text}\n\n{prompts[analysis_type]}"}
                    ],
                    temperature=0.3,
                    max_tokens=2000
                )

                logger.debug("OpenAI response received successfully")
                print("OpenAI response received successfully")  # Debug print

                if not response or not response.choices:
                    error_msg = "Resposta vazia da OpenAI"
                    logger.error(error_msg)
                    print(error_msg)  # Debug print
                    st.error(error_msg)
                    return None

                analysis_result = response.choices[0].message.content
                if not analysis_result:
                    error_msg = "Análise retornou vazia"
                    logger.error(error_msg)
                    print(error_msg)  # Debug print
                    st.error(error_msg)
                    return None

                logger.debug(f"Analysis result length: {len(analysis_result)}")
                print(f"Analysis result length: {len(analysis_result)}")  # Debug print
                return analysis_result

            except Exception as api_error:
                error_msg = f"Erro na API da OpenAI: {str(api_error)}"
                logger.error(error_msg, exc_info=True)
                print(error_msg)  # Debug print
                st.error(error_msg)
                return None

        except Exception as e:
            error_msg = f"Erro durante a análise: {str(e)}"
            logger.error(error_msg, exc_info=True)
            print(error_msg)  # Debug print
            st.error(error_msg)
            return None


def get_api_keys():
    """Obtém as chaves API das variáveis de ambiente"""
    openai_key = os.getenv("OPENAI_API_KEY")
    pubmed_key = os.getenv("NCBI_API_KEY")

    if not openai_key and 'OPENAI_API_KEY' in st.secrets:
        openai_key = st.secrets['OPENAI_API_KEY']
    if not pubmed_key and 'NCBI_API_KEY' in st.secrets:
        pubmed_key = st.secrets['NCBI_API_KEY']

    return openai_key, pubmed_key


def main():
    st.set_page_config(
        page_title="PubMed Searcher",
        page_icon="🔍",
        layout="wide"
    )

    st.title("🔍 Buscador e Analisador PubMed")
    st.write("Pesquise e analise artigos científicos do PubMed")

    # Initialize session state if not exists
    if 'analysis_results' not in st.session_state:
        st.session_state.analysis_results = {}

    openai_key, pubmed_key = get_api_keys()

    with st.sidebar:
        st.header("Configurações")
        email = st.text_input("Email", key="email")

        pubmed_key_input = st.text_input(
            "API Token PubMed",
            value=pubmed_key if pubmed_key else "",
            type="password",
            key="pubmed_key"
        )
        openai_key_input = st.text_input(
            "API Token OpenAI",
            value=openai_key if openai_key else "",
            type="password",
            key="openai_key"
        )

        pubmed_key = pubmed_key_input or pubmed_key
        openai_key = openai_key_input or openai_key

        st.markdown("---")
        st.markdown("""
        ### Dicas de Pesquisa
        - Use AND para combinar termos
        - Use OR para alternativas
        - Use aspas para frases exatas
        - Use [MeSH Terms] para termos MeSH
        """)

    with st.form("search_form"):
        col1, col2 = st.columns([3, 1])
        with col1:
            query = st.text_input("Termo de busca",
                                  placeholder='Ex: "diabetes mellitus"[MeSH Terms]')
        with col2:
            max_results = st.number_input("Máx. resultados",
                                          min_value=1, max_value=100, value=10)

        col3, col4 = st.columns(2)
        with col3:
            start_year = st.number_input("Ano inicial",
                                         min_value=1900,
                                         max_value=datetime.now().year,
                                         value=2020)
        with col4:
            end_year = st.number_input("Ano final",
                                       min_value=1900,
                                       max_value=datetime.now().year,
                                       value=datetime.now().year)

        submitted = st.form_submit_button("Buscar")

    if submitted and email and pubmed_key and query:
        logger.debug("Search form submitted")
        searcher = PubMedSearch(email, pubmed_key, openai_key)
        results = searcher.search_articles(query, max_results, start_year, end_year)

        # Modifique a parte do código que processa os resultados da busca:

    if results:
        logger.debug(f"Found {len(results)} results")
        df = pd.DataFrame(results)

        # Criar uma seção para mostrar análises dos artigos gratuitos
        st.markdown("## 📊 Análise dos Artigos Gratuitos")

        # Filtrar apenas artigos gratuitos
        free_articles = df[df['is_free'] == True]
        total_articles = len(df)
        free_count = len(free_articles)

        st.info(f"Encontrados {free_count} artigos gratuitos de um total de {total_articles} artigos.")

        if free_count > 0:
            with st.spinner("Analisando artigos gratuitos..."):
                analyses = []

                # Progress bar
                progress_bar = st.progress(0)

                for idx, article in free_articles.iterrows():
                    try:
                        # Atualizar progresso
                        progress = (idx + 1) / free_count
                        progress_bar.progress(progress)

                        st.text(f"Processando artigo {idx + 1} de {free_count}: {article['title']}")

                        # Obter texto completo
                        full_text = searcher.get_full_article(article['pmc_id'])

                        if full_text:
                            # Realizar análise
                            analysis = searcher.analyze_article(full_text, "complete")

                            if analysis:
                                analyses.append({
                                    'title': article['title'],
                                    'authors': article['authors'],
                                    'year': article['year'],
                                    'pmid': article['pmid'],
                                    'analysis': analysis
                                })

                    except Exception as e:
                        logger.error(f"Erro ao processar artigo {article['pmid']}: {str(e)}")
                        continue

                # Remover barra de progresso
                progress_bar.empty()

                # Mostrar resultados
                if analyses:
                    st.success(f"✅ Análise completa de {len(analyses)} artigos gratuitos!")

                    # Criar tabs para cada artigo analisado
                    tabs = st.tabs([f"Artigo {i+1}" for i in range(len(analyses))])

                    for idx, (tab, analysis) in enumerate(zip(tabs, analyses)):
                        with tab:
                            st.markdown(f"### {analysis['title']}")
                            st.markdown(f"**Autores:** {analysis['authors']}")
                            st.markdown(f"**Ano:** {analysis['year']} | **PMID:** {analysis['pmid']}")

                            with st.expander("Ver Análise Completa", expanded=True):
                                st.markdown(analysis['analysis'])

                            # Adicionar botão de download para cada análise
                            st.download_button(
                                "📥 Download desta análise",
                                analysis['analysis'],
                                f"analise_{analysis['pmid']}.txt",
                                "text/plain",
                                key=f"download_{analysis['pmid']}"
                            )

                    # Botão para download de todas as análises em um único arquivo
                    all_analyses = "\n\n" + "="*50 + "\n\n".join(
                        f"TÍTULO: {a['title']}\n"
                        f"AUTORES: {a['authors']}\n"
                        f"ANO: {a['year']} | PMID: {a['pmid']}\n\n"
                        f"ANÁLISE:\n{a['analysis']}"
                        for a in analyses
                    )

                    st.download_button(
                        "📥 Download de todas as análises",
                        all_analyses,
                        f"todas_analises_{int(time.time())}.txt",
                        "text/plain",
                        key='download-all-analyses'
                    )

                    # Gerar resumo geral
                    st.markdown("## 📋 Resumo Geral")
                    summary_prompt = f"""
                    Com base nas análises dos {len(analyses)} artigos acima, forneça:
                    1. Principais temas e descobertas comuns
                    2. Tendências metodológicas
                    3. Lacunas e limitações comuns
                    4. Implicações práticas gerais
                    5. Recomendações para pesquisas futuras
                    """

                    try:
                        all_analyses_text = "\n\n".join(a['analysis'] for a in analyses)
                        general_summary = searcher.analyze_article(all_analyses_text, "complete")
                        if general_summary:
                            st.markdown(general_summary)

                            st.download_button(
                                "📥 Download do resumo geral",
                                general_summary,
                                f"resumo_geral_{int(time.time())}.txt",
                                "text/plain",
                                key='download-summary'
                            )
                    except Exception as e:
                        st.error(f"Erro ao gerar resumo geral: {str(e)}")

                else:
                    st.warning("Não foi possível analisar nenhum dos artigos gratuitos.")

        # Mostrar todos os resultados originais
        st.markdown("## 📚 Todos os Resultados")
        for idx, article in df.iterrows():
            with st.container():
                st.markdown(f"""
                ### {article['title']}
                **Autores:** {article['authors']}  
                **Journal:** {article['journal']} ({article['year']})  
                **Abstract:** {article['abstract']}
                
                [Ver no PubMed]({article['url']}) | [DOI](https://doi.org/{article['doi']})
                """)

                if article['is_free']:
                    st.success("✅ Artigo gratuito disponível!")

                st.markdown("---")

        st.download_button(
            "📥 Download todos os resultados (CSV)",
            df.to_csv(index=False).encode('utf-8'),
            f"pubmed_results_{int(time.time())}.csv",
            "text/csv",
            key='download-csv'
        )

    else:
        st.warning("Nenhum resultado encontrado.")
        logger.warning("No results found for the search")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error("Application error:", exc_info=True)
        st.error(f"Erro na aplicação: {str(e)}")
