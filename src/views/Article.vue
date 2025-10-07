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
        <div class="card shadow-sm border-0 bg-white rounded-3">
          <div class="card-body p-4">
            <div class="article-content">
          <div v-if="currentPost" class="article-meta pb-2 mb-0">
            <h1 class="article-title mb-3">{{ currentPost.title }}</h1>
            <div class="text-secondary d-flex flex-wrap gap-3 align-items-center">
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
    <!-- 移动端：TOC抽屉按钮和右侧抽屉（全局一致样式，从右侧弹出） -->
    <TocDrawerButton v-if="rawMarkdown" />
    <MobileTocDrawer v-if="rawMarkdown" />
  </div>
</template>

<script>
import RenderMarkdown from '@/components/layout/RenderMarkdown.vue'
import OnThisPage from '@/components/layout/OnThisPage.vue'
import TocDrawerButton from '@/components/widgets/TocDrawerButton.vue'
import MobileTocDrawer from '@/components/widgets/MobileTocDrawer.vue'
import NavigationTree from '@/components/layout/NavigationTree.vue'
import notesFlat from '@/content/notes.json'
import projectsFlat from '@/content/projects.json'
import topicsFlat from '@/content/topics.json'

const markdownModules = import.meta.glob('../content-src/**/*.md', { query: '?raw', import: 'default', eager: false });
const assetModules = import.meta.glob('../content-src/**/*.{png,jpg,jpeg,gif,svg,webp}', { as: 'url', eager: true });
const keys = Object.keys(markdownModules);

export default {
  name: 'ArticleView',
  components: { RenderMarkdown, OnThisPage, NavigationTree, TocDrawerButton, MobileTocDrawer },
  props: { path: { type: [String, Array], default: '' } },
  data() {
    return { rawMarkdown: '', showLeft: false, showRight: false, currentPath: '', allArticles: [], groupedArticles: {} }
  },
  computed: {
    isDesktop() { return window.innerWidth >= 992 },
    isProject() {
      return !!(this.currentPost?.path && this.currentPost.path.startsWith('projects/'));
    },
    isNote() {
      return !!(this.currentPost?.path && this.currentPost.path.startsWith('notes/'));
    },
    currentPost() {
      if (!this.currentPath) return null;
      return this.allArticles.find(article => article.path.replace(/\.md$/, '') === this.currentPath.replace(/\.md$/, ''));
    },
    currentCategoryArticles() {
      return this.currentPost ? this.groupedArticles[this.currentPost.category] || [] : [];
    },
    currentCategoryIndex() {
      return this.currentPost ? this.currentCategoryArticles.findIndex(article => article.path === this.currentPost.path) : -1;
    },
    prevPost() {
      return this.currentCategoryIndex > 0 ? this.currentCategoryArticles[this.currentCategoryIndex - 1] : null;
    },
    nextPost() {
      const lastIdx = this.currentCategoryArticles.length - 1;
      return this.currentCategoryIndex >= 0 && this.currentCategoryIndex < lastIdx ? this.currentCategoryArticles[this.currentCategoryIndex + 1] : null;
    },
    readingMinutes() {
      if (!this.rawMarkdown) return 0;
      const text = this.rawMarkdown.replace(/(```[\s\S]*?```)|(!\[.*?\]\(.*?\))|(\[.*?\]\(.*?\))|([#>*\-]+)/g, '').trim();
      return Math.max(1, Math.round(text.length / 800));
    }
  },
  created() {
    this.initAllArticles();
    this.loadArticleContent();
  },
  mounted() {
    this.showLeft = this.showRight = this.isDesktop;
    window.addEventListener('resize', this.onResize);
    window.addEventListener('scroll', this.updateSidebarDimensions);
    window.addEventListener('resize', this.updateSidebarDimensions);
  },
  watch: {
    '$route': {
      handler(to, from) {
        const oldPath = this.normalizeRoutePathParam(from?.params?.path);
        const newPath = this.normalizeRoutePathParam(to?.params?.path);
        if (oldPath !== newPath) {
          this.$refs.onThisPageRef?.resetToc();
          this.loadArticleContent();
          // 移动端：路由变更后自动关闭左侧抽屉
          if (this.showMobileSidebar) this.closeMobileSidebar();
        }
      },
      immediate: false, deep: true
    },
    rawMarkdown() {
      this.$nextTick(() => this.updateSidebarDimensions());
    }
  },
  beforeUnmount() {
    window.removeEventListener('resize', this.onResize);
    window.removeEventListener('scroll', this.updateSidebarDimensions);
    window.removeEventListener('resize', this.updateSidebarDimensions);
  },
  methods: {
    toArticle(p) {
      return { name: 'Article', params: { path: p.replace(/\.md$/, '').split('/') } };
    },

    // 移动端与桌面端不同处理：移动端打开统一抽屉，桌面端切换本页侧栏
    onMobileOrDesktopToggle(side) {
      if (!this.isDesktop) {
        // 通知 Header 打开移动端抽屉（带灰色半透明遮罩）
        window.dispatchEvent(new Event('open-mobile-sidebar'));
        return;
      }
      if (side === 'left') this.showLeft = !this.showLeft;
      if (side === 'right') this.showRight = !this.showRight;
    },
    iconFor(side) {
      if (!this.isDesktop) {
        // 移动端始终显示打开图标，抽屉内使用叉关闭
        return side === 'left' ? 'bi bi-list' : 'bi bi-bookmark';
      }
      return side === 'left'
        ? (this.showLeft ? 'bi bi-x-lg' : 'bi bi-list')
        : (this.showRight ? 'bi bi-x-lg' : 'bi bi-bookmark');
    },
    labelFor(side) {
      if (!this.isDesktop) {
        return side === 'left' ? '显示目录' : '显示本页目录';
      }
      return side === 'left'
        ? (this.showLeft ? '隐藏目录' : '显示目录')
        : (this.showRight ? '隐藏本页目录' : '显示本页目录');
    },

    onResize() {
      this.showLeft = this.showRight = this.isDesktop;
      this.updateSidebarDimensions();
    },
    getMatchedKey(rel) {
      // 文章路径为 notes/xxx/xxx.md，真实文件在 content-src/notes 下
      const normalized = rel.replace(/^notes\//, 'notes/');
      const candidates = [
        `/content-src/${normalized}`,
        `${normalized}`,
        `/content-src/${normalized}?raw`,
        `${normalized}?raw`,
        `../content-src/${normalized}`,
        `../content-src/${normalized}?raw`
      ];
      return keys.find(k => candidates.some(suf => k.endsWith(suf)));
    },
    normalizeRoutePathParam(p) {
      if (Array.isArray(p)) return p.join('/');
      if (typeof p === 'string' && p.length) return p;
      return '';
    },
    initAllArticles() {
      const allArticlesList = [];
      const groupedArticlesMap = {};

      // 名称映射（与 contentGenerator 保持一致）
      const nameMap = {
        'Programming': '编程语言',
        'Bioinformatics': '生物信息学',
        'Omics': '组学技术',
        'DataScience': '数据科学',
        'python': 'Python',
        'r': 'R语言',
        'shell': 'Shell',
        'javascript': 'JavaScript',
        'alignment': '序列比对',
        'structure': '结构分析',
        'genomics': '基因组学',
        'proteomics': '蛋白质组学',
        'transcriptomics': '转录组学',
        'statistics': '统计分析',
        'machinelearning': '机器学习',
        'visualization': '数据可视化'
      };
      const formatName = (n) => nameMap[n] || n;

      // 从扁平 notes.json 生成 Article 使用的条目
      if (Array.isArray(notesFlat)) {
        notesFlat.forEach(item => {
          const pathWithMd = `notes/${item.relativePath}.md`;
          const seg0 = item.relativePath.split('/')[0] || 'notes';
          const categoryName = formatName(seg0);
          const art = {
            title: item.title || item.relativePath.split('/').pop(),
            path: pathWithMd,
            date: item.date || '',
            tags: item.tags || [],
            preview: item.description || '',
            category: categoryName
          };
          if (!groupedArticlesMap[categoryName]) groupedArticlesMap[categoryName] = [];
          groupedArticlesMap[categoryName].push(art);
          allArticlesList.push(art);
        });
      }

      // 从 projects.json 生成项目条目（仅标题，日期/标签为空）
      if (Array.isArray(projectsFlat)) {
        projectsFlat.forEach(item => {
          const pathWithMd = `projects/${item.relativePath}.md`;
          const seg0 = item.relativePath.split('/')[0] || 'projects';
          const categoryName = formatName(seg0);
          const art = {
            title: item.title || item.relativePath.split('/').pop(),
            path: pathWithMd,
            date: '',
            tags: [],
            preview: '',
            category: categoryName
          };
          if (!groupedArticlesMap[categoryName]) groupedArticlesMap[categoryName] = [];
          groupedArticlesMap[categoryName].push(art);
          allArticlesList.push(art);
        });
      }

      // 从 topics.json 生成主题条目（同项目：无日期/标签，仅标题）
      if (Array.isArray(topicsFlat)) {
        topicsFlat.forEach(item => {
          const pathWithMd = `topics/${item.relativePath}.md`;
          const seg0 = item.relativePath.split('/')[0] || 'topics';
          const categoryName = formatName(seg0);
          const art = {
            title: item.title || item.relativePath.split('/').pop(),
            path: pathWithMd,
            date: '',
            tags: [],
            preview: '',
            category: categoryName
          };
          if (!groupedArticlesMap[categoryName]) groupedArticlesMap[categoryName] = [];
          groupedArticlesMap[categoryName].push(art);
          allArticlesList.push(art);
        });
      }

      this.allArticles = allArticlesList;
      this.groupedArticles = groupedArticlesMap;
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
        const rewritten = this.rewriteImageLinks(markdownModule, this.currentPath);
        this.rawMarkdown = rewritten;

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
      const [header, footer, leftContent, rightContent, leftContainer, rightContainer] = [
        document.querySelector('header'),
        document.querySelector('footer'),
        this.$refs.leftSidebarContent,
        this.$refs.rightSidebarContent,
        this.$refs.leftSidebarContainer,
        this.$refs.rightSidebarContainer
      ];
      if (!leftContent || !rightContent || !leftContainer || !rightContainer) return;

      const [headerH, footerH, viewportH, scrollTop, docH] = [
        header?.offsetHeight || 0,
        footer?.offsetHeight || 0,
        window.innerHeight,
        window.scrollY,
        document.documentElement.scrollHeight
      ];
      const remainingH = Math.max(0, docH - scrollTop - headerH - footerH - 40);
      const availableH = Math.min(viewportH - headerH - 40, remainingH);

      [leftContent, rightContent].forEach(el => {
        el.style.maxHeight = `${availableH}px`;
        el.style.overflowY = el.scrollHeight > availableH ? 'auto' : 'visible';
      });
    },

    // 将 Markdown 中的相对图片路径重写为可访问的打包 URL
    // articlePath 形如 'notes/Programming/javascript/vue/vue.md'
    rewriteImageLinks(md, articlePath) {
      try {
        const articleDir = articlePath.replace(/^[./]*/, '').replace(/\.md$/, '').split('/').slice(0, -1).join('/');
        // 在 assetModules 中，key 形如 '../content-src/notes/.../image.png'
        const toAssetUrl = (relPath) => {
          // 忽略外链与以 / 开头的绝对路径
          if (/^(https?:)?\/\//i.test(relPath) || relPath.startsWith('/')) return relPath;
          // 归一化 ./ ../
          const parts = (articleDir + '/' + relPath).split('/').filter(p => p && p !== '.');
          const stack = [];
          parts.forEach(p => {
            if (p === '..') stack.pop();
            else stack.push(p);
          });
          const normalized = stack.join('/');
          const candidateKey = `../content-src/${normalized}`;
          const url = assetModules[candidateKey];
          return url || relPath; // 找不到则保留原样
        };

        // 替换图片语法与 HTML img
        const mdReplaced = md
          .replace(/!\[([^\]]*)\]\(([^)]+)\)/g, (m, alt, src) => {
            const clean = src.trim().replace(/^<|>|&/g, '');
            const url = toAssetUrl(clean);
            return `![${alt}](${url})`;
          })
          .replace(/<img\s+([^>]*?)src=["']([^"']+)["'](.*?)>/gi, (m, pre, src, post) => {
            const url = toAssetUrl(src.trim());
            return `<img ${pre}src="${url}"${post}>`;
          });

        return mdReplaced;
      } catch (e) {
        console.warn('rewriteImageLinks failed', e);
        return md;
      }
    },
    handleMarkdownRendered() {
      this.enhanceHeadings();
      this.updateSidebarDimensions();
      this.$nextTick(() => this.$refs.onThisPageRef?.refreshToc());
    },
    enhanceHeadings() {
      const container = this.$refs.mainContent?.querySelector('.markdown-body');
      if (!container) return;
      this.cleanDuplicateH1(container);
      this.addAnchorLinks(container);
    },
    cleanDuplicateH1(container) {
      if (!this.currentPost) return;
      const pageTitle = this.currentPost.title.trim().toLowerCase();
      container.querySelectorAll('h1').forEach(h1 => {
        const h1Text = h1.textContent.trim().toLowerCase();
        if (h1Text === pageTitle) h1.remove();
        else {
          const h2 = document.createElement('h2');
          Array.from(h1.attributes).forEach(attr => h2.setAttribute(attr.name, attr.value));
          h2.innerHTML = h1.innerHTML;
          h1.replaceWith(h2);
        }
      });
    },
    addAnchorLinks(container) {
      const scrollToHeading = (heading) => {
        const targetTop = window.scrollY + heading.getBoundingClientRect().top - 8;
        window.scrollTo({ top: Math.max(0, targetTop), behavior: 'smooth' });
      };

      container.querySelectorAll('h2, h3, h4, h5, h6').forEach(heading => {
        heading.querySelector('.heading-anchor')?.remove();
        const anchorBtn = document.createElement('button');
        anchorBtn.type = 'button';
        anchorBtn.className = 'heading-anchor';
        anchorBtn.setAttribute('aria-label', '置顶当前标题');
        anchorBtn.setAttribute('tabindex', '0');
        anchorBtn.setAttribute('aria-hidden', 'false');
        anchorBtn.textContent = '#';

        anchorBtn.addEventListener('click', (e) => { e.stopPropagation(); scrollToHeading(heading); });
        anchorBtn.addEventListener('keydown', (e) => {
          if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); scrollToHeading(heading); }
        });

        heading.appendChild(anchorBtn);
      });
    }
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
.article-meta { padding-bottom: 0.75rem !important; }

.article-title { font-size: 2.1rem; font-weight: 700; }
.article-meta .text-secondary { font-size: 1rem; }
.article-meta .badge { font-size: 0.95rem; font-weight: 500; }
.article-meta .tag-badge { 
  color: #212529 !important; 
  background-color: #f1f3f5 !important; 
}

:deep(.markdown-body) {
  font-size: 1.125rem; line-height: 1.8; color: var(--bs-gray-800);
}
:deep(.markdown-body p) { margin-bottom: 0.75rem; }

:deep(.markdown-body h2, .markdown-body h3, .markdown-body h4, .markdown-body h5, .markdown-body h6) {
  position: relative; padding-right: 0; display: inline-block; width: 100%;
  margin-top: 0.5rem; margin-bottom: 0.5rem;
  padding-bottom: 0.25rem;
}
:deep(.markdown-body h2) {
  font-size: 1.8rem;
  font-weight: 700;
  margin-top: 0.5rem;
  border-bottom: 1px solid var(--bs-border-color);
}
:deep(.markdown-body h3) {
  font-size: 1.6rem;
  font-weight: 700;
}
:deep(.markdown-body h4) {
  font-size: 1.4rem;
  font-weight: 600;
}
:deep(.markdown-body h5) {
  font-size: 1.25rem;
  font-weight: 600;
}
:deep(.markdown-body h6) {
  font-size: 1.125rem;
  font-weight: 600;
}

:deep(.heading-anchor) {
  position: relative;
  display: inline-block;
  margin-left: 0.3em;
  color: var(--bs-gray-600);
  font-size: 0.9em;
  font-weight: 400;
  opacity: 0.6;
  transition: opacity 0.2s ease, color 0.2s ease;
  cursor: pointer;
  background: none;
  border: none;
  border-radius: 0;
  box-shadow: none;
  padding: 0;
  line-height: 1;
}
:deep(
  .markdown-body h2:hover .heading-anchor,
  .markdown-body h3:hover .heading-anchor,
  .markdown-body h4:hover .heading-anchor,
  .markdown-body h5:hover .heading-anchor,
  .markdown-body h6:hover .heading-anchor,
  .heading-anchor:focus
) {
  opacity: 1;
  color: var(--bs-primary);
}
:deep(.heading-anchor:focus) {
  outline: none;
  outline-offset: 0;
}

.article-navigation {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 0.75rem;
  margin-top: 2rem;
  border-top: 1px solid var(--bs-border-color);
  padding-top: 1rem;
}
.article-nav-item {
  display: flex;
  align-items: center;
  gap: 1rem;
  border: 1px solid var(--bs-border-color);
  border-radius: .5rem;
  padding: 1rem;
  text-decoration: none;
  color: var(--bs-body-color);
  transition: all 0.2s ease-in-out;
  min-width: 0;
}
.article-nav-item:hover {
  border-color: var(--bs-primary);
  color: var(--bs-primary);
  transform: translateY(-2px);
  box-shadow: 0 4px 12px rgba(var(--bs-primary-rgb), 0.1);
}
.article-nav-item.prev { grid-column: 1; }
.article-nav-item.next { grid-column: 2; justify-content: flex-end; text-align: right; }
.nav-arrow {
  font-size: 1.5rem; font-weight: 300; line-height: 1;
  color: var(--bs-gray-400); transition: color 0.2s ease-in-out;
}
.article-nav-item:hover .nav-arrow { color: var(--bs-primary); }
.nav-details { display: flex; flex-direction: column; overflow: hidden; min-width: 0; }

.article-nav-item.next .nav-details { align-items: flex-end; flex: 1 1 auto; }
.nav-label {
  font-size: 0.875rem; color: var(--bs-gray-600);
  transition: color 0.2s ease-in-out;
}
.article-nav-item:hover .nav-label { color: var(--bs-primary); }
.nav-title { font-weight: 500; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; direction: ltr; }
.article-nav-item.next .nav-title { text-align: left; max-width: 100%; }

@media (max-width: 991px) {
  .sticky-sidebar {
    position: static; top: auto; bottom: auto !important;
    max-height: none !important; overflow-y: visible !important;
  }
  .navigation-container { margin-bottom: 1rem; }
  .toc-container { margin-top: 1rem; }
  .article-content { padding: 1rem; }
}


/* 移动端侧栏按钮：与 Header 扁平风格一致 */
.sidebar-toggle {
  background: transparent !important;
  border: none !important;
  color: var(--bs-body-color);
  transition: color 0.15s ease-in-out, transform 0.1s ease-in-out;
}
.sidebar-toggle:hover,
.sidebar-toggle:focus {
  background-color: transparent !important;
  color: var(--bs-primary);
}
.sidebar-toggle i {
  font-size: 1.05rem;
  vertical-align: middle;
}

/* 保持与 Header 一致的交互反馈 */
.sidebar-toggle:active {
  transform: scale(0.98);
}

</style>
