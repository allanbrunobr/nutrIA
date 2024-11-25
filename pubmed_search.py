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

class PubMedSearch:
    def __init__(self, email: str, pubmed_key: str, openai_key: str = None):
        self.email = email
        self.pubmed_key = pubmed_key
        self.openai_key = openai_key
        Entrez.email = email
        Entrez.api_key = pubmed_key
        if openai_key:
            self.openai_client = OpenAI(api_key=openai_key)

    def search_articles(self,
                        query: str,
                        max_results: int = 10,
                        start_year: int = None,
                        end_year: int = None) -> List[Dict]:
        """Realiza busca no PubMed"""
        if start_year and end_year:
            query = f"{query} AND ({start_year}:{end_year}[pdat])"

        try:
            with st.spinner('Buscando artigos...'):
                search_handle = Entrez.esearch(
                    db="pubmed",
                    term=query,
                    retmax=max_results,
                    sort="relevance"
                )
                search_results = Entrez.read(search_handle)
                search_handle.close()
                logging.debug("DEBUG: Iniciando análise do artigo")
                if not search_results["IdList"]:
                    return []

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

                        processed_articles.append({
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
                        })

                    except Exception as e:
                        st.error(f"Erro ao processar artigo: {str(e)}")
                        continue

                return processed_articles

        except Exception as e:
            st.error(f"Erro na busca: {str(e)}")
            return []

    def get_full_article(self, pmc_id: str) -> str:
        """Obtém o texto completo do artigo do PMC"""
        try:
            handle = Entrez.efetch(
                db="pmc",
                id=pmc_id,
                rettype="xml",
                retmode="xml"
            )
            xml_content = handle.read()
            root = ET.fromstring(xml_content)

            article_text = []

            title = root.find(".//article-title")
            if title is not None:
                article_text.append(f"TÍTULO: {title.text}")

            abstract = root.find(".//abstract")
            if abstract is not None:
                article_text.append("\nRESUMO:")
                for p in abstract.findall(".//p"):
                    if p.text:
                        article_text.append(p.text)

            body = root.find(".//body")
            if body is not None:
                article_text.append("\nCONTEÚDO PRINCIPAL:")
                for section in body.findall(".//sec"):
                    title = section.find("title")
                    if title is not None:
                        article_text.append(f"\n{title.text}:")
                    for p in section.findall(".//p"):
                        if p.text:
                            article_text.append(p.text)

            return "\n".join(article_text)
        except Exception as e:
            st.error(f"Erro ao extrair conteúdo do PMC: {str(e)}")
            return ""

    def analyze_article(self, text: str, analysis_type: str = "complete") -> str:
        """Analisa o artigo usando OpenAI"""
        logging.debug("DEBUG: Iniciando análise do artigo")

        if not self.openai_key:
            return "Chave da OpenAI não configurada."

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
            print("DEBUG: Resposta do OpenAI recebida")
            return response.choices[0].message.content
        except Exception as e:
            st.error(f"Erro na análise: {str(e)}")
            return ""


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
    st.title("🔍 Buscador e Analisador PubMed")
    st.write("Pesquise e analise artigos científicos do PubMed")

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
        searcher = PubMedSearch(email, pubmed_key, openai_key)
        results = searcher.search_articles(query, max_results, start_year, end_year)

        if results:
            df = pd.DataFrame(results)

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

                        article_key = f"analyze_{article['pmid']}"

                        if article_key not in st.session_state:
                            st.session_state[article_key] = {
                                "full_text": None,
                                "analysis": None,
                                "show_analysis": False
                            }

                        if st.button(f"Carregar artigo {article['pmid']}", key=f"load_{article_key}"):
                            with st.spinner("Obtendo texto completo do artigo..."):
                                full_text = searcher.get_full_article(article['pmc_id'])
                                print(f"DEBUG: Texto completo carregado: {bool(full_text)}")
                                if full_text:
                                    st.session_state[article_key]["full_text"] = full_text
                                    st.session_state[article_key]["show_analysis"] = True
                                else:
                                    st.error("Não foi possível obter o texto completo do artigo.")

                        if st.session_state[article_key]["show_analysis"]:
                            full_text = st.session_state[article_key]["full_text"]
                            if full_text:
                                st.markdown("### Texto Completo do Artigo")
                                st.text(full_text[:1000])

                                analysis_type = st.radio(
                                    "Tipo de análise:",
                                    ["complete", "methodology", "results", "practical"],
                                    key=f"type_{article_key}"
                                )

                                if st.button("Executar Análise", key=f"execute_analysis_{article_key}"):
                                    print("DEBUG: Botão Executar Análise clicado")
                                    with st.spinner("Analisando artigo..."):
                                        analysis = searcher.analyze_article(full_text, analysis_type)
                                        print(f"DEBUG: Análise realizada? {'Sim' if analysis else 'Não'}")
                                        if analysis:
                                            st.session_state[article_key]["analysis"] = analysis
                                            st.markdown("### Resultado da Análise")
                                            st.markdown(analysis)
                                            st.download_button(
                                                "📥 Download da análise",
                                                analysis,
                                                f"analise_{article['pmid']}.txt",
                                                "text/plain",
                                                key=f"download_analysis_{article_key}"
                                            )
                                        else:
                                            st.error("Erro durante a análise.")

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

    elif submitted:
        st.error("Por favor, preencha email, API token e termo de busca.")


if __name__ == "__main__":
    main()
