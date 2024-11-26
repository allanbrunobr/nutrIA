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

if 'search_results' not in st.session_state:
    st.session_state.search_results = None
if 'articles_text' not in st.session_state:
    st.session_state.articles_text = {}
if 'searcher' not in st.session_state:
    st.session_state.searcher = None

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

    # Na parte do analyze_article, modifique o dicionário de prompts para corresponder aos tipos de análise:

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
            logger.debug("Preparing OpenAI request")

            # Mapeia os tipos de análise para prompts
            prompts = {
                "Análise Completa": """
                    Por favor, forneça uma análise completa deste artigo científico, incluindo:
                    1. Principais descobertas e conclusões
                    2. Metodologia utilizada
                    3. Limitações do estudo
                    4. Implicações práticas
                    5. Recomendações baseadas nas evidências
                """,
                "Resumo Principal": "Resuma os pontos principais e descobertas mais importantes deste artigo científico.",
                "Metodologia": "Analise detalhadamente a metodologia usada neste estudo, destacando pontos fortes e limitações.",
                "Resultados e Conclusões": "Apresente uma análise detalhada dos resultados e conclusões do estudo.",
                "Implicações Práticas": """
                    Foque nas implicações práticas deste estudo:
                    1. Como os resultados podem ser aplicados na prática
                    2. Recomendações práticas
                    3. Pontos de atenção para implementação
                """,
                "Crítica do Estudo": """
                    Faça uma análise crítica do estudo, considerando:
                    1. Qualidade metodológica
                    2. Validade dos resultados
                    3. Potenciais vieses
                    4. Lacunas na pesquisa
                    5. Sugestões para estudos futuros
                """
            }

            if analysis_type not in prompts:
                logger.error(f"Tipo de análise inválido: {analysis_type}")
                print(f"Tipo de análise inválido: {analysis_type}")  # Debug print
                return None

            # Limitar o tamanho do texto
            max_text_length = 14000
            text = text[:max_text_length] + "..." if len(text) > max_text_length else text

            try:
                logger.debug(f"Enviando requisição para análise do tipo: {analysis_type}")
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

                if not response or not response.choices:
                    error_msg = "Resposta vazia da OpenAI"
                    logger.error(error_msg)
                    st.error(error_msg)
                    return None

                analysis_result = response.choices[0].message.content
                logger.debug(f"Análise concluída com sucesso: {len(analysis_result)} caracteres")
                return analysis_result

            except Exception as api_error:
                error_msg = f"Erro na API da OpenAI: {str(api_error)}"
                logger.error(error_msg, exc_info=True)
                st.error(error_msg)
                return None

        except Exception as e:
            error_msg = f"Erro durante a análise: {str(e)}"
            logger.error(error_msg, exc_info=True)
            st.error(error_msg)
            return None

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

                        # Verifica se é artigo gratuito e extrai PMCID
                        is_free = False
                        pmcid = None
                        if 'ArticleIdList' in article['PubmedData']:
                            for id_item in article['PubmedData']['ArticleIdList']:
                                if hasattr(id_item, 'attributes') and id_item.attributes.get('IdType') == 'pmc':
                                    is_free = True
                                    pmcid = str(id_item)
                                    break

                        # Extrai DOI
                        doi = None
                        if 'ELocationID' in article_data:
                            for id in article_data['ELocationID']:
                                if hasattr(id, 'attributes') and id.attributes.get('EIdType') == 'doi':
                                    doi = str(id)
                                    break

                        # Se não encontrou DOI nos ELocationID, procura no ArticleIdList
                        if not doi and 'ArticleIdList' in article['PubmedData']:
                            for id_item in article['PubmedData']['ArticleIdList']:
                                if hasattr(id_item, 'attributes') and id_item.attributes.get('IdType') == 'doi':
                                    doi = str(id_item)
                                    break

                        pmid = str(article['MedlineCitation']['PMID'])

                        processed_article = {
                            'title': article_data['ArticleTitle'],
                            'authors': '; '.join(authors[:3]) + (' et al.' if len(authors) > 3 else ''),
                            'journal': article_data['Journal']['Title'],
                            'year': int(article_data['Journal']['JournalIssue']['PubDate'].get('Year', 0)),
                            'abstract': abstract[:300] + '...' if len(abstract) > 300 else abstract,
                            'pmid': pmid,
                            'pmcid': pmcid,  # Agora incluímos o PMCID
                            'doi': doi if doi else 'N/A',
                            'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                            'is_free': is_free,
                            'pmc_url': f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/" if pmcid else None
                        }

                        processed_articles.append(processed_article)
                        logger.debug(f"Processed article PMID: {pmid} | PMCID: {pmcid if pmcid else 'None'}")

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

    def get_full_article(self, pmcid: str) -> str:
        """Obtém o texto completo do artigo do PMC usando o PMCID"""
        logger.debug(f"Iniciando get_full_article para PMCID: {pmcid}")
        print(f"Iniciando get_full_article para PMCID: {pmcid}")  # Debug print

        if not pmcid:
            error_msg = "PMCID está vazio ou None"
            logger.error(error_msg)
            st.error(error_msg)
            return ""

        try:
            # Remove prefixo PMC se existir
            pmcid = pmcid.replace("PMC", "")
            logger.debug(f"PMCID processado: {pmcid}")

            try:
                import requests
                from bs4 import BeautifulSoup

                # URL correta do PMC
                url = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{pmcid}"
                logger.debug(f"Tentando acessar URL: {url}")
                print(f"Tentando acessar URL: {url}")  # Debug print

                headers = {
                    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36'
                }

                response = requests.get(url, headers=headers)

                if response.status_code == 200:
                    soup = BeautifulSoup(response.text, 'html.parser')
                    article_text = []

                    # Tenta encontrar o título
                    title = soup.find('h1', {'class': 'content-title'})
                    if title:
                        article_text.append(f"TÍTULO: {title.text.strip()}")

                    # Tenta encontrar o abstract
                    abstract = soup.find('div', {'class': 'abstract'})
                    if abstract:
                        article_text.append("\nRESUMO:")
                        abstract_text = abstract.get_text().strip()
                        article_text.append(abstract_text)

                    # Tenta encontrar o corpo do artigo
                    article_body = soup.find('div', {'class': 'article-body'})
                    if article_body:
                        article_text.append("\nCONTEÚDO PRINCIPAL:")

                        # Processa seções
                        sections = article_body.find_all(['h2', 'h3', 'p'])
                        for section in sections:
                            if section.name in ['h2', 'h3']:
                                article_text.append(f"\n{section.text.strip()}:")
                            else:
                                text = section.get_text().strip()
                                if text:
                                    article_text.append(text)

                    # Se não encontrou o corpo do artigo, tenta outra abordagem
                    if not article_body:
                        main_content = soup.find('article') or soup.find('div', {'class': 'main-content'})
                        if main_content:
                            for p in main_content.find_all('p'):
                                text = p.get_text().strip()
                                if text:
                                    article_text.append(text)

                    if article_text:
                        full_text = "\n\n".join(article_text)
                        logger.debug(f"Texto extraído com sucesso. Tamanho: {len(full_text)}")
                        return full_text
                    else:
                        error_msg = "Não foi possível extrair o texto do artigo"
                        logger.error(error_msg)
                        st.error(error_msg)
                        return ""
                else:
                    error_msg = f"Erro ao acessar o artigo. Status code: {response.status_code}"
                    logger.error(error_msg)
                    st.error(error_msg)
                    return ""

            except Exception as web_error:
                error_msg = f"Erro ao acessar página do PMC: {str(web_error)}"
                logger.error(error_msg)
                st.error(error_msg)
                return ""

        except Exception as e:
            error_msg = f"Erro geral: {str(e)}"
            logger.error(error_msg)
            st.error(error_msg)
            return ""

    # Modificação na parte que processa os resultados da busca:
def process_article_metadata(self, article_data: dict) -> dict:
    """Processa os metadados do artigo"""
    try:
        pmid = str(article_data['MedlineCitation']['PMID'])
        logger.debug(f"Processando metadados para PMID: {pmid}")
        # Processa autores
        # Processa autores
        authors = []
        if 'AuthorList' in article_data['MedlineCitation']['Article']:
            for author in article_data['MedlineCitation']['Article']['AuthorList']:
                if 'LastName' in author and 'ForeName' in author:
                    authors.append(f"{author['LastName']} {author['ForeName']}")

        # Processa abstract
        abstract = ""
        if 'Abstract' in article_data['MedlineCitation']['Article']:
            abstract_text = article_data['MedlineCitation']['Article']['Abstract']['AbstractText']
            if isinstance(abstract_text, list):
                abstract = ' '.join(str(text) for text in abstract_text)
            else:
                abstract = str(abstract_text)

        # Processa identificadores
        pmcid = None
        doi = "N/A"
        is_free = False

        # Verifica DOI
        if 'ELocationID' in article_data['MedlineCitation']['Article']:
            for id in article_data['MedlineCitation']['Article']['ELocationID']:
                if hasattr(id, 'attributes') and id.attributes.get('EIdType') == 'doi':
                    doi = str(id)
                    break

        # Verifica PMC ID e status free
        if 'ArticleIdList' in article_data.get('PubmedData', {}):
            for id_item in article_data['PubmedData']['ArticleIdList']:
                if hasattr(id_item, 'attributes'):
                    id_type = id_item.attributes.get('IdType', '')
                    if id_type == 'pmc':
                        pmcid = str(id_item)
                        is_free = True
                        break
                    elif id_type == 'doi' and not doi:
                        doi = str(id_item)

        # Cria URL do PMC se disponível
        pmc_url = None
        if pmcid:
            clean_pmcid = pmcid.replace('PMC', '')
            pmc_url = f"https://pmc.ncbi.nlm.nih.gov/articles/PMC{clean_pmcid}"
            logger.debug(f"URL do PMC gerada: {pmc_url}")

        return {
            'title': article_data['ArticleTitle'],
            'authors': '; '.join(authors[:3]) + (' et al.' if len(authors) > 3 else ''),
            'journal': article_data['Journal']['Title'],
            'year': int(article_data['Journal']['JournalIssue']['PubDate'].get('Year', 0)),
            'abstract': abstract[:300] + '...' if len(abstract) > 300 else abstract,
            'pmid': pmid,
            'pmcid': pmcid,
            'doi': doi if doi else 'N/A',
            'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            'is_free': is_free,
            'pmc_url': pmc_url
        }

    except Exception as e:
        logger.error(f"Erro ao processar metadados do artigo: {str(e)}")
        raise




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
        st.session_state.searcher = PubMedSearch(email, pubmed_key, openai_key)
        results = st.session_state.searcher.search_articles(query, max_results, start_year, end_year)
        st.session_state.search_results = results

    if st.session_state.search_results:
        results = st.session_state.search_results
        df = pd.DataFrame(results)


        if results:
            logger.debug(f"Found {len(results)} results")
            df = pd.DataFrame(results)

            # Mostrar resumo dos resultados
            st.markdown("## 📊 Resumo dos Resultados")
            total_articles = len(df)
            free_articles = df[df['is_free'] == True]
            free_count = len(free_articles)

            col1, col2 = st.columns(2)
            with col1:
                st.metric("Total de Artigos", total_articles)
            with col2:
                st.metric("Artigos Gratuitos", free_count)

            # Processar artigos gratuitos
            if free_count > 0:
                st.markdown("## 📚 Artigos Gratuitos")

                for idx, article in free_articles.iterrows():
                    with st.expander(f"📄 {article['title']}", expanded=False):
                        st.markdown(f"""
                        **Autores:** {article['authors']}  
                        **Journal:** {article['journal']} ({article['year']})  
                        **Abstract:** {article['abstract']}
                        
                        [Ver no PubMed]({article['url']}) | [DOI](https://doi.org/{article['doi']})
                        """)

                        if article['pmc_url'] and article['is_free']:
                            st.success("✅ Artigo gratuito disponível!")
                            st.markdown(f"[Ver artigo completo no PMC]({article['pmc_url']})")

                        # Botões para carregar e analisar
                        col1, col2 = st.columns(2)
                        with col1:
                            article_key = f"article_{article['pmid']}"

                            # Botão para mostrar URL
                            if st.button("🔗 Ver URL do Artigo", key=f"url_{article['pmid']}"):
                                pmc_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{article['pmcid']}"
                                st.info(f"""
                                **URL do artigo:**  
                                {pmc_url}
                                """)
                                st.code(pmc_url, language=None)
                                st.markdown(f"[Abrir artigo no PMC]({pmc_url})")

                            # Botão para carregar texto
                            if st.button("📥 Carregar Texto Completo", key=f"load_{article['pmid']}"):
                                try:
                                    # Recria o searcher se necessário
                                    if st.session_state.searcher is None:
                                        st.session_state.searcher = PubMedSearch(email, pubmed_key, openai_key)

                                    with st.spinner("Obtendo texto completo..."):
                                        logger.debug(f"Tentando carregar artigo {article['pmcid']}")
                                        pmc_url = f"https://pmc.ncbi.nlm.nih.gov/articles/{article['pmcid']}"
                                        st.info(f"Acessando: {pmc_url}")

                                        full_text = st.session_state.searcher.get_full_article(article['pmcid'])

                                        if full_text:
                                            # Salva o texto na session_state
                                            st.session_state.articles_text[article_key] = full_text
                                            st.success("✅ Texto carregado com sucesso!")

                                            # Mostrar prévia do texto
                                            st.markdown("### Prévia do Artigo")
                                            st.text(full_text[:500] + "...")

                                        else:
                                            st.error("❌ Não foi possível carregar o texto completo.")
                                            logger.error(f"Texto vazio retornado para {article['pmcid']}")

                                except Exception as e:
                                    error_msg = f"Erro ao carregar artigo: {str(e)}"
                                    logger.error(error_msg)
                                    st.error(f"❌ {error_msg}")

                        with col2:
                            # Verifica se o texto está carregado
                            if article_key in st.session_state.articles_text:

                                if st.session_state.searcher is None:
                                    st.session_state.searcher = PubMedSearch(email, pubmed_key, openai_key)

                                # Adiciona select box para tipo de análise
                                analysis_type = st.selectbox(
                                    "Tipo de Análise:",
                                    [
                                        "Análise Completa",
                                        "Resumo Principal",
                                        "Metodologia",
                                        "Resultados e Conclusões",
                                        "Implicações Práticas",
                                        "Crítica do Estudo"
                                    ],
                                    key=f"analysis_type_{article['pmid']}"
                                )

                                # Mapeia os tipos de análise para prompts específicos
                                analysis_prompts = {
                                    "Análise Completa": """
                                        Por favor, forneça uma análise completa deste artigo científico, incluindo:
                                        1. Principais descobertas e conclusões
                                        2. Metodologia utilizada
                                        3. Limitações do estudo
                                        4. Implicações práticas
                                        5. Recomendações baseadas nas evidências
                                    """,
                                    "Resumo Principal": "Resuma os pontos principais e descobertas mais importantes deste artigo científico.",
                                    "Metodologia": "Analise detalhadamente a metodologia usada neste estudo, destacando pontos fortes e limitações.",
                                    "Resultados e Conclusões": "Apresente uma análise detalhada dos resultados e conclusões do estudo.",
                                    "Implicações Práticas": """
                                        Foque nas implicações práticas deste estudo:
                                        1. Como os resultados podem ser aplicados na prática
                                        2. Recomendações práticas
                                        3. Pontos de atenção para implementação
                                    """,
                                    "Crítica do Estudo": """
                                        Faça uma análise crítica do estudo, considerando:
                                        1. Qualidade metodológica
                                        2. Validade dos resultados
                                        3. Potenciais vieses
                                        4. Lacunas na pesquisa
                                        5. Sugestões para estudos futuros
                                    """
                                }

                                # Mostrar prévia do texto carregado (sem usar expander)
                                st.markdown("### Texto Carregado (Prévia)")
                                st.text(st.session_state.articles_text[article_key][:500] + "...")

                                # Botão para análise
                                if st.button("🤖 Analisar com ChatGPT", key=f"analyze_{article['pmid']}"):
                                    with st.spinner(f"Realizando {analysis_type}..."):
                                        try:
                                            analysis = st.session_state.searcher.analyze_article(
                                                st.session_state.articles_text[article_key],
                                                analysis_type
                                            )
                                            if analysis:
                                                st.success("✅ Análise concluída!")
                                                st.markdown(f"### Análise do Artigo ({analysis_type})")
                                                st.markdown(analysis)

                                                # Botões para download
                                                col3, col4 = st.columns([1, 1])
                                                with col3:
                                                    st.download_button(
                                                        "📥 Download da análise",
                                                        analysis,
                                                        f"analise_{article['pmid']}_{analysis_type.lower().replace(' ', '_')}.txt",
                                                        "text/plain",
                                                        key=f"download_{article['pmid']}_{analysis_type}"
                                                    )
                                                with col4:
                                                    st.download_button(
                                                        "📥 Download do texto original",
                                                        st.session_state.articles_text[article_key],
                                                        f"texto_original_{article['pmid']}.txt",
                                                        "text/plain",
                                                        key=f"download_original_{article['pmid']}"
                                                    )
                                            else:
                                                st.error("❌ Erro ao gerar análise.")
                                        except Exception as e:
                                            st.error(f"❌ Erro durante análise: {str(e)}")
                            else:
                                st.info("Carregue o texto completo do artigo para realizar a análise.")

            # Mostrar todos os resultados
            st.markdown("## 📑 Todos os Resultados")
            for _, article in df.iterrows():
                with st.expander(f"{article['title']}", expanded=False):
                    st.markdown(f"""
                    **Autores:** {article['authors']}  
                    **Journal:** {article['journal']} ({article['year']})  
                    **Abstract:** {article['abstract']}
                    
                    [Ver no PubMed]({article['url']}) | [DOI](https://doi.org/{article['doi']})
                    """)

                    if article['is_free']:
                        st.success("✅ Artigo gratuito disponível!")
                        if article['pmcid']:
                            st.markdown(f"[Ver artigo completo no PMC](https://www.ncbi.nlm.nih.gov/pmc/articles/{article['pmcid']})")

            # Botão para download de todos os resultados
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
    else:
        if submitted:
            st.error("Por favor, preencha email, API token e termo de busca.")
            logger.error("Form submitted with missing required fields")

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logger.error("Application error:", exc_info=True)
        st.error(f"Erro na aplicação: {str(e)}")
