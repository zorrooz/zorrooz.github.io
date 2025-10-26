<template>
  <div class="container">
    <div class="row py-4 px-0">
      <!-- 左侧：根目录（2列） -->
      <div class="col-12 col-lg-2 order-2 order-lg-1" ref="leftSidebarContainer">
        <div class="sticky-sidebar" ref="leftSidebarContent">
          <div v-if="isDesktop" class="navigation-container">
            <NavigationTree />
          </div>
        </div>
      </div>

      <!-- 中间：正文（8列） -->
      <div class="col-12 col-lg-8 order-1 order-lg-2" ref="mainContent">
        <div class="card shadow-sm border-0 rounded-3" :style="{ backgroundColor: 'var(--app-card-bg)' }">
          <div class="card-body p-4">
            <div class="article-content">
              <div v-if="currentPost" class="article-meta pb-2 mb-0">
                <h1 class="article-title mb-3">{{ currentPost.title }}</h1>
                <div class="d-flex flex-wrap gap-3 align-items-center" :style="{ color: 'var(--app-text-muted)' }">
                  <span v-if="isNote && currentPost.date"><i class="bi bi-calendar3 me-1"></i>更新于 {{ currentPost.date }}</span>
                  <span v-if="readingMinutes"><i class="bi bi-clock me-1"></i>阅读约 {{ readingMinutes }} 分钟</span>
                </div>
                <div v-if="currentPost.tags?.length" class="d-flex flex-wrap gap-2 mt-2">
                  <span v-for="(tag, idx) in currentPost.tags" :key="idx" class="badge tag-badge fw-normal py-1 px-2 rounded-3">
                    # {{ tag }}
                  </span>
                </div>
              </div>

              <RenderMarkdown 
                v-if="rawMarkdown" 
                :rawMarkdown="rawMarkdown" 
                :articlePath="currentPost?.path || ''"
                :articleTitle="currentPost?.title || ''"
                @markdown-rendered="handleMarkdownRendered"
              />

              <!-- 上一篇 / 下一篇 -->
              <nav class="article-navigation" v-if="rawMarkdown">
                <router-link v-if="prevPost" :to="toArticle(prevPost.path)" class="article-nav-item prev">
                  <div class="nav-arrow"><</div>
                  <div class="nav-details">
                    <div class="nav-label">上一页</div>
                    <div class="nav-title">{{ prevPost.title }}</div>
                  </div>
                </router-link>
                <router-link v-if="nextPost" :to="toArticle(nextPost.path)" class="article-nav-item next">
                  <div class="nav-details">
                    <div class="nav-label">下一页</div>
                    <div class="nav-title">{{ nextPost.title }}</div>
                  </div>
                  <div class="nav-arrow">></div>
                </router-link>
              </nav>
            </div>
          </div>
        </div>
      </div>

      <!-- 右侧：本页目录（2列） -->
      <div class="col-12 col-lg-2 order-3" ref="rightSidebarContainer">
        <div class="sticky-sidebar" ref="rightSidebarContent">
          <div v-if="isDesktop" class="toc-container">
            <OnThisPage 
              ref="onThisPageRef"
              containerSelector=".markdown-body" 
              :levels="[2, 3]" 
            />
          </div>
        </div>
      </div>
    </div>

    <!-- 移动端：TOC抽屉 -->
    <TocDrawer v-if="rawMarkdown" />
  </div>
</template>

<script>
import RenderMarkdown from '@/components/layout/RenderMarkdown.vue'
import OnThisPage from '@/components/layout/OnThisPage.vue'
import TocDrawer from '@/components/widgets/TocDrawer.vue'
import NavigationTree from '@/components/layout/NavigationTree.vue'
import categoryData from '@/content/categories.json'

const markdownModules = import.meta.glob('../content-src/**/*.md', { query: '?raw', import: 'default', eager: false });
const keys = Object.keys(markdownModules);

export default {
  name: 'ArticleView',
  components: { RenderMarkdown, OnThisPage, NavigationTree, TocDrawer },
  props: { path: { type: [String, Array], default: '' } },
  data() {
    return { 
      rawMarkdown: '', 
      currentPath: '', 
      allArticles: [], 
      groupedArticles: {}, 
      viewportWidth: (typeof window !== 'undefined' ? window.innerWidth : 1024) 
    }
  },
  computed: {
    isDesktop() { return this.viewportWidth >= 992 },
    isNote() { return !!(this.currentPost?.path && this.currentPost.path.startsWith('notes/')); },
    
    currentPost() {
      if (!this.currentPath) return null;
      return this.allArticles.find(article => article.path.replace(/\.md$/, '') === this.currentPath.replace(/\.md$/, ''));
    },

    groupLinearArticles() {
      if (!this.currentPost) return [];
      const [type, group] = this.currentPost.path.replace(/\.md$/, '').split('/');
      const linear = [];
      
      const pushFromUrl = (title, articleUrl) => {
        if (!articleUrl) return;
        const parts = String(articleUrl).replace(/^\/+/, '').split('/');
        const i0 = parts[0] === 'article' ? 1 : 0;
        const t = parts[i0], g = parts[i0 + 1];
        if (t !== type || g !== group) return;
        const rest = parts.slice(i0 + 2);
        const pathNoExt = `${t}/${g}/${rest.join('/')}`;
        linear.push({ title, path: `${pathNoExt}.md` });
      };

      if (Array.isArray(categoryData)) {
        for (const section of categoryData) {
          if (!Array.isArray(section.items)) continue;
          for (const item of section.items) {
            if (item?.name !== group) continue;
            if (Array.isArray(item.articles)) {
              item.articles.forEach(a => pushFromUrl(a.title, a.articleUrl));
            }
            if (Array.isArray(item.categories)) {
              item.categories.forEach(cat => {
                if (Array.isArray(cat.articles)) {
                  cat.articles.forEach(a => pushFromUrl(a.title, a.articleUrl));
                }
              });
            }
          }
        }
      }
      return linear;
    },

    currentLinearIndex() {
      if (!this.currentPost) return -1;
      return this.groupLinearArticles.findIndex(a => a.path.replace(/\.md$/, '') === this.currentPost.path.replace(/\.md$/, ''));
    },

    prevPost() {
      const idx = this.currentLinearIndex;
      return idx > 0 ? this.groupLinearArticles[idx - 1] : null;
    },

    nextPost() {
      const idx = this.currentLinearIndex;
      const last = this.groupLinearArticles.length - 1;
      return idx >= 0 && idx < last ? this.groupLinearArticles[idx + 1] : null;
    },

    readingMinutes() {
      if (!this.rawMarkdown) return 0;
      const text = this.rawMarkdown.replace(/(```[\s\S]*?```)|(!\[.*?\]\(.*?\))|(\[.*?\]\(.*?\))|([#>*\-]+)/g, '').trim();
      return Math.max(1, Math.round(text.length / 800));
    }
  },
  created() {
    this.buildFromCategories();
    this.loadArticleContent();
  },
  mounted() {
    this.viewportWidth = window.innerWidth;
    window.addEventListener('resize', this.onResize);
    window.addEventListener('scroll', this.updateSidebarDimensions);
  },
  watch: {
    '$route': {
      handler(to, from) {
        const oldPath = this.normalizeRoutePathParam(from?.params?.path);
        const newPath = this.normalizeRoutePathParam(to?.params?.path);
        if (oldPath !== newPath) {
          this.$refs.onThisPageRef?.resetToc();
          this.loadArticleContent();
        }
      },
      immediate: false
    },
    rawMarkdown() {
      this.$nextTick(() => this.updateSidebarDimensions());
    }
  },
  beforeUnmount() {
    window.removeEventListener('resize', this.onResize);
    window.removeEventListener('scroll', this.updateSidebarDimensions);
  },
  methods: {
    toArticle(p) {
      return { name: 'Article', params: { path: p.replace(/\.md$/, '').split('/') } };
    },

    buildFromCategories() {
      const all = [], grouped = {};

      const pushArticle = (artTitle, articleUrl, tags = [], dateStr = '') => {
        if (typeof articleUrl !== 'string' || !articleUrl.trim()) return;
        const parts = articleUrl.replace(/^\/+/, '').split('/');
        const idxArticle = parts[0] === 'article' ? 1 : 0;
        const t = parts[idxArticle], g = parts[idxArticle + 1];
        const restParts = parts.slice(idxArticle + 2);
        const subKey = restParts[0] || '__root__';
        const rest = restParts.join('/');
        const pathNoExt = `${t}/${g}/${rest}`;

        const art = { title: artTitle, path: `${pathNoExt}.md`, date: dateStr || '', tags: Array.isArray(tags) ? tags : [], preview: '', category: `${t}/${g}/${subKey}` };
        if (!grouped[art.category]) grouped[art.category] = [];
        grouped[art.category].push(art);
        all.push(art);
      };

      if (Array.isArray(categoryData)) {
        categoryData.forEach(section => {
          if (!Array.isArray(section.items)) return;
          section.items.forEach(item => {
            const itemLatest = item?.stats?.latestDate || '';
            if (Array.isArray(item.articles)) {
              item.articles.forEach(a => pushArticle(a.title, a.articleUrl, a?.tags || [], itemLatest));
            }
            if (Array.isArray(item.categories)) {
              item.categories.forEach(cat => {
                const catLatest = cat?.stats?.latestDate || itemLatest || '';
                if (Array.isArray(cat.articles)) {
                  cat.articles.forEach(a => pushArticle(a.title, a.articleUrl, a?.tags || [], catLatest));
                }
              });
            }
          });
        });
      }

      this.allArticles = all;
      this.groupedArticles = grouped;
    },

    onResize() {
      this.viewportWidth = window.innerWidth;
      this.updateSidebarDimensions();
    },

    getMatchedKey(rel) {
      const normalized = rel.replace(/^notes\//, 'notes/');
      const candidates = [
        `/content-src/${normalized}`, `${normalized}`,
        `/content-src/${normalized}?raw`, `${normalized}?raw`,
        `../content-src/${normalized}`, `../content-src/${normalized}?raw`
      ];
      return keys.find(k => candidates.some(suf => k.endsWith(suf)));
    },

    normalizeRoutePathParam(p) {
      return Array.isArray(p) ? p.join('/') : (typeof p === 'string' && p.length ? p : '');
    },

    async loadArticleContent() {
      this.rawMarkdown = '';
      try {
        const currentPathClean = this.normalizeRoutePathParam(this.$route.params.path);
        const matchedPost = this.allArticles.find(article => article.path.replace(/\.md$/, '') === currentPathClean);
        if (!matchedPost) throw new Error(`Article not found: ${currentPathClean}`);
        this.currentPath = matchedPost.path;

        const matchedKey = this.getMatchedKey(this.currentPath);
        if (!matchedKey) throw new Error(`Markdown not found: ${this.currentPath}`);

        const markdownModule = await markdownModules[matchedKey]();
        this.rawMarkdown = markdownModule;

        this.$nextTick(() => {
          this.updateSidebarDimensions();
          this.$refs.onThisPageRef?.refreshToc();
        });
      } catch (error) {
        console.error('Load article failed:', error);
        this.rawMarkdown = '# Article Not Found\n\nThe requested article could not be loaded. Please check the URL.';
        this.$nextTick(() => this.$refs.onThisPageRef?.refreshToc());
      }
    },

    updateSidebarDimensions() {
      const header = document.querySelector('header'), footer = document.querySelector('footer');
      const leftContent = this.$refs.leftSidebarContent, rightContent = this.$refs.rightSidebarContent;
      if (!leftContent || !rightContent || !this.$refs.leftSidebarContainer || !this.$refs.rightSidebarContainer) return;

      const headerH = header?.offsetHeight || 0, footerH = footer?.offsetHeight || 0;
      const viewportH = window.innerHeight, scrollTop = window.scrollY, docH = document.documentElement.scrollHeight;
      const remainingH = Math.max(0, docH - scrollTop - headerH - footerH - 40);
      const availableH = Math.min(viewportH - headerH - 40, remainingH);

      [leftContent, rightContent].forEach(el => {
        el.style.maxHeight = `${availableH}px`;
        el.style.overflowY = el.scrollHeight > availableH ? 'auto' : 'visible';
      });
    },

    handleMarkdownRendered() {
      this.updateSidebarDimensions();
      this.$nextTick(() => this.$refs.onThisPageRef?.refreshToc());
    },
  }
}
</script>

<style scoped>
.sticky-sidebar {
  position: sticky; top: 30px; box-sizing: border-box; width: 100%;
  -webkit-overflow-scrolling: touch; transition: max-height 0.2s ease;
}
.navigation-container { margin-bottom: 0; }
.toc-container { margin-top: 0; }
.article-content { min-height: 400px; padding: 0; }
.article-meta { padding-bottom: 0.75rem !important; border-bottom: none !important; }

.article-title { font-size: 2.1rem; font-weight: 700; color: var(--app-text-emphasis); }
.article-meta .text-secondary { font-size: 1rem; color: var(--app-text-muted); }
.article-meta .badge { font-size: 0.95rem; font-weight: 500; }
.article-meta .tag-badge { 
  color: var(--app-tag-text) !important; 
  background-color: var(--app-tag-bg) !important; 
}

:deep(.markdown-body) {
  font-size: 1.125rem; line-height: 1.8; color: var(--app-text);
}
:deep(.markdown-body p) { margin-bottom: 0.75rem; }



.article-navigation {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 0.75rem;
  margin-top: 2rem;
  padding-top: 1rem;
  border-top: none !important;
}
.article-nav-item {
  display: flex;
  align-items: center;
  gap: 1rem;
  border: 1px solid var(--app-border);
  border-radius: .5rem;
  padding: 1rem;
  text-decoration: none;
  color: var(--app-text);
  transition: all 0.2s ease-in-out;
  min-width: 0;
}
.article-nav-item:hover {
  border-color: var(--app-primary);
  color: var(--app-primary);
  transform: translateY(-2px);
  box-shadow: 0 4px 12px var(--app-primary-rgb-01);
}
.article-nav-item.prev { grid-column: 1; }
.article-nav-item.next { grid-column: 2; justify-content: flex-end; text-align: right; }
.nav-arrow { font-size: 1.5rem; font-weight: 300; line-height: 1; color: var(--app-nav-arrow-color); transition: color 0.2s ease-in-out; }
.article-nav-item:hover .nav-arrow { color: var(--app-primary); }
.nav-details { display: flex; flex-direction: column; overflow: hidden; min-width: 0; }
.article-nav-item.next .nav-details { align-items: flex-end; flex: 1 1 auto; }
.nav-label { font-size: 0.875rem; color: var(--app-text-muted); transition: color 0.2s ease-in-out; }
.article-nav-item:hover .nav-label { color: var(--app-primary); }
.nav-title { font-weight: 500; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; direction: ltr; }
.article-nav-item.next .nav-title { text-align: left; max-width: 100%; }

@media (max-width: 991px) {
  .sticky-sidebar { position: static; top: auto; bottom: auto !important; max-height: none !important; overflow-y: visible !important; }
  .navigation-container { margin-bottom: 1rem; }
  .toc-container { margin-top: 1rem; }
  .card-body { padding: 0.75rem !important; }
  .article-content { padding: 0.25rem; }
}
</style>
