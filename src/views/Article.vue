<template>
  <div class="container">
    <div class="row pt-4 px-0">
      <!-- 左侧：根目录（2列） -->
      <div class="col-12 col-lg-2 order-2 order-lg-1 px-2" ref="leftSidebarContainer">
        <div class="sticky-sidebar" ref="leftSidebarContent">
          <button class="btn btn-sm btn-outline-secondary d-lg-none mb-2 w-100" @click="showLeft = !showLeft">
            <i class="bi bi-list"></i> {{ showLeft ? '隐藏目录' : '显示目录' }}
          </button>
          <div v-show="showLeft || isDesktop" class="navigation-container">
            <NavigationTree />
          </div>
        </div>
      </div>

      <!-- 中间：正文（8列） -->
      <div class="col-12 col-lg-8 order-1 order-lg-2 px-0" ref="mainContent">
        <div class="article-content">
          <div v-if="currentPost" class="article-meta border-bottom pb-3 mb-3">
            <h1 class="h3 mb-2">{{ currentPost.title }}</h1>
            <div class="text-secondary small d-flex flex-wrap gap-3 align-items-center">
              <span><i class="bi bi-calendar3 me-1"></i>更新于 {{ currentPost.date }}</span>
              <span v-if="readingMinutes"><i class="bi bi-clock me-1"></i>阅读约 {{ readingMinutes }} 分钟</span>
            </div>
            <div v-if="currentPost.tags?.length" class="d-flex flex-wrap gap-2 mt-2">
              <span v-for="(tag, idx) in currentPost.tags" :key="idx" class="badge bg-white text-body fw-normal py-1 px-2 small rounded-3">
                # {{ tag }}
              </span>
            </div>
          </div>

          <RenderMarkdown :rawMarkdown="rawMarkdown" />

          <!-- 上一篇 / 下一篇（优化后） -->
          <nav class="article-nav border-top pt-4 mt-5 d-flex flex-wrap justify-content-between align-items-center gap-3">
            <router-link
              v-if="prevPost"
              :to="toArticle(prevPost.path)"
              class="article-nav-btn btn btn-light border rounded-3 px-4 py-2 text-body d-flex align-items-center gap-2 transition-all duration-200"
              :aria-label="`上一篇：${prevPost.title}`"
            >
              <i class="bi bi-arrow-left"></i>
              <span class="article-nav-text">上一篇</span>
              <span class="article-nav-title d-none d-sm-inline-block">：{{ prevPost.title }}</span>
            </router-link>
            <span v-else class="article-nav-empty text-secondary small italic px-4 py-2">没有更早的文章</span>

            <router-link
              v-if="nextPost"
              :to="toArticle(nextPost.path)"
              class="article-nav-btn btn btn-light border rounded-3 px-4 py-2 text-body d-flex align-items-center gap-2 transition-all duration-200 ms-auto"
              :aria-label="`下一篇：${nextPost.title}`"
            >
              <span class="article-nav-text">下一篇</span>
              <span class="article-nav-title d-none d-sm-inline-block">：{{ nextPost.title }}</span>
              <i class="bi bi-arrow-right"></i>
            </router-link>
            <span v-else class="article-nav-empty text-secondary small italic px-4 py-2 ms-auto">没有更新的文章</span>
          </nav>
        </div>
      </div>

      <!-- 右侧：本页目录（2列） -->
      <div class="col-12 col-lg-2 order-3 px-2" ref="rightSidebarContainer">
        <div class="sticky-sidebar" ref="rightSidebarContent">
          <button class="btn btn-sm btn-outline-secondary d-lg-none mb-2 w-100" @click="showRight = !showRight">
            <i class="bi bi-bookmark"></i> {{ showRight ? '隐藏本页目录' : '显示本页目录' }}
          </button>
          <div v-show="showRight || isDesktop" class="toc-container">
            <OnThisPage containerSelector=".markdown-body" :levels="[2, 3]" />
          </div>
        </div>
      </div>
    </div>
  </div>
</template>

<script>
import RenderMarkdown from '@/components/layout/RenderMarkdown.vue'
import OnThisPage from '@/components/layout/OnThisPage.vue'
import NavigationTree from '@/components/layout/NavigationTree.vue'
import posts from '@/content/posts.json'

const markdownModules = import.meta.glob('../content/**/*.md', { as: 'raw' });
const keys = Object.keys(markdownModules);

export default {
  name: 'ArticleView',
  components: { RenderMarkdown, OnThisPage, NavigationTree },
  props: { path: { type: [String, Array], default: '' } },
  data() {
    return { rawMarkdown: '', showLeft: false, showRight: false, currentPath: '', currentIndex: -1 }
  },
  computed: {
    isDesktop() { return window.innerWidth >= 992 },
    currentPost() { return this.currentIndex >= 0 ? posts[this.currentIndex] : null },
    prevPost() { return this.currentIndex - 1 >= 0 ? posts[this.currentIndex - 1] : null },
    nextPost() { return this.currentIndex + 1 < posts.length ? posts[this.currentIndex + 1] : null },
    readingMinutes() {
      if (!this.rawMarkdown) return 0
      const text = this.rawMarkdown.replace(
        /(```[\s\S]*?```)|(!\[.*?\]\(.*?\))|(\[.*?\]\(.*?\))|([#>*\-]+)/g, 
        ''
      ).trim()
      return Math.max(1, Math.round(text.length / 800))
    }
  },
  mounted() {
    this.showLeft = this.isDesktop
    this.showRight = this.isDesktop
    window.addEventListener('resize', this.onResize)
    window.addEventListener('scroll', this.updateSidebarDimensions)
    window.addEventListener('resize', this.updateSidebarDimensions)
    
    this.loadArticleContent()
    this.$nextTick(() => this.updateSidebarDimensions())
  },
  watch: {
    '$route': {
      handler(to, from) {
        const oldPath = from?.params?.path?.join('/')
        const newPath = to?.params?.path?.join('/')
        if (oldPath !== newPath) this.loadArticleContent()
      },
      immediate: true, deep: true
    }
  },
  beforeUnmount() {
    window.removeEventListener('resize', this.onResize)
    window.removeEventListener('scroll', this.updateSidebarDimensions)
    window.removeEventListener('resize', this.updateSidebarDimensions)
  },
  methods: {
    toArticle(p) { return { name: 'Article', params: { path: p } } },
    onResize() {
      this.showLeft = this.isDesktop
      this.showRight = this.isDesktop
      this.updateSidebarDimensions()
    },
    getMatchedKey(rel) {
      const suffixes = [
        `/content/${rel}`, rel, `/content/${rel}?raw`, `${rel}?raw`,
        `../content/${rel}`, `../content/${rel}?raw`
      ]
      return keys.find(k => suffixes.some(suf => k.endsWith(suf)))
    },
    async loadArticleContent() {
      console.log('Available markdown modules:', keys);
      try {
        const routeParam = this.$route.params.path
        const markdownPath = Array.isArray(routeParam) 
          ? routeParam.join('/') 
          : (routeParam || this.path || 'notes/Programming/python/25-09-19--python-fastq.md')
        const targetRel = markdownPath.endsWith('.md') ? markdownPath : `${markdownPath}.md`
        console.log('Attempting to load article:', targetRel);

        const matchedKey = this.getMatchedKey(targetRel)
        if (!matchedKey) throw new Error(`No matching module for: ${targetRel}`)

        const mod = await markdownModules[matchedKey]()
        this.rawMarkdown = typeof mod === 'string' ? mod : (mod.default ?? mod)
        this.currentPath = targetRel
        this.currentIndex = posts.findIndex(p => p.path === this.currentPath)
        this.$nextTick(() => setTimeout(this.enhanceHeadings, 0))
      } catch (error) {
        console.error('Load failed:', error);
        const fallbackRel = 'notes/Programming/python/25-09-19--python-fastq.md'
        const fallbackKey = this.getMatchedKey(fallbackRel)
        try {
          if (!fallbackKey) throw new Error('Fallback key not found')
          const mod = await markdownModules[fallbackKey]()
          this.rawMarkdown = typeof mod === 'string' ? mod : (mod.default ?? mod)
          this.currentPath = fallbackRel
          this.currentIndex = posts.findIndex(p => p.path === this.currentPath)
          this.$nextTick(() => setTimeout(this.enhanceHeadings, 0))
        } catch (fallbackError) {
          console.error('Fallback load failed:', fallbackError);
          this.rawMarkdown = '# Article Loading Failed\n\nPlease check the article path.'
        }
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
      ]
      if (!leftContent || !rightContent || !leftContainer || !rightContainer) return

      const [headerH, footerH, viewportH, scrollTop, docH] = [
        header?.offsetHeight || 0,
        footer?.offsetHeight || 0,
        window.innerHeight,
        window.scrollY,
        document.documentElement.scrollHeight
      ]
      const remainingH = Math.max(0, docH - scrollTop - headerH - footerH - 40)
      const availableH = Math.min(viewportH - headerH - 40, remainingH)

      [leftContent, rightContent].forEach(el => {
        el.style.maxHeight = `${availableH}px`
        el.style.overflowY = el.scrollHeight > availableH ? 'auto' : 'visible'
      })
    },
    enhanceHeadings() {
      const container = this.$refs.mainContent?.querySelector('.markdown-body')
      if (!container) return

      Array.from(container.querySelectorAll('h1')).forEach(h1 => {
        const h2 = document.createElement('h2')
        Array.from(h1.attributes).forEach(attr => h2.setAttribute(attr.name, attr.value))
        h2.innerHTML = h1.innerHTML
        h1.replaceWith(h2)
      })

      try {
        const topTitle = this.currentPost?.title?.trim()
        const firstHeading = container.querySelector('h1, h2, h3, h4, h5, h6')
        if (topTitle && firstHeading?.textContent?.trim() === topTitle) firstHeading.remove()
      } catch (e) {}

      const scrollToTop = (el) => {
        const targetTop = Math.max(0, window.scrollY + el.getBoundingClientRect().top - 8)
        window.scrollTo({ top: targetTop, behavior: 'smooth' })
      }

      Array.from(container.querySelectorAll('h2, h3, h4, h5, h6')).forEach(h => {
        if (!h.querySelector('.heading-anchor')) {
          const btn = Object.assign(document.createElement('button'), {
            type: 'button',
            className: 'heading-anchor',
            'aria-label': '置顶',
            tabindex: '0',
            'aria-hidden': 'true',
            textContent: ''
          })
          btn.addEventListener('click', (e) => { e.stopPropagation(); scrollToTop(h) })
          btn.addEventListener('keydown', (e) => {
            if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); scrollToTop(h) }
          })
          h.appendChild(btn)
        }
      })
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
.article-content { min-height: 400px; padding: 0 1.5rem; }

/* 标题与锚点样式 */
:deep(.markdown-body h2, .markdown-body h3, .markdown-body h4, .markdown-body h5, .markdown-body h6) {
  position: relative; padding-right: 0; display: inline-block; width: 100%;
}
:deep(.heading-anchor) {
  position: relative; display: inline-block; margin-left: 0.3em;
  color: var(--bs-gray-400); font-size: 0.9em; font-weight: 400;
  opacity: 0; transition: opacity 0.2s ease; cursor: pointer;
  background: none; border: none; border-radius: 0; box-shadow: none;
}
:deep(.heading-anchor::before) { content: '#'; }
:deep(
  .markdown-body h2:hover .heading-anchor,
  .markdown-body h3:hover .heading-anchor,
  .markdown-body h4:hover .heading-anchor,
  .markdown-body h5:hover .heading-anchor,
  .markdown-body h6:hover .heading-anchor,
  .heading-anchor:focus
) { opacity: 1; }
:deep(.heading-anchor:hover) { color: var(--bs-gray-500); opacity: 1; }
:deep(.heading-anchor:focus) { outline: 1px solid transparent; outline-offset: 2px; }

/* 上一篇/下一篇导航样式 */
.article-nav {
  background-color: var(--bs-gray-50);
  border-radius: 8px;
  padding: 1rem;
  margin-top: 2rem !important;
}

.article-nav-btn {
  &:hover {
    background-color: var(--bs-primary-50);
    border-color: var(--bs-primary-200);
    color: var(--bs-primary);
    transform: translateY(-1px);
  }
  &:active {
    background-color: var(--bs-primary-100);
    border-color: var(--bs-primary-300);
    transform: translateY(0);
  }
  &:focus-visible {
    box-shadow: 0 0 0 3px rgba(var(--bs-primary-rgb), 0.2);
    border-color: var(--bs-primary);
    outline: none;
  }
}

.article-nav-title {
  max-width: 180px;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
}

.article-nav-empty {
  height: 100%;
  display: flex;
  align-items: center;
  color: var(--bs-gray-500) !important;
}

/* 响应式样式 */
@media (max-width: 991px) {
  .sticky-sidebar {
    position: static; top: auto; bottom: auto !important;
    max-height: none !important; overflow-y: visible !important;
  }
  .navigation-container { margin-bottom: 1rem; }
  .toc-container { margin-top: 1rem; }
  .article-content { padding: 0 0.5rem; }
  :deep(.markdown-body h2, .markdown-body h3, .markdown-body h4, .markdown-body h5, .markdown-body h6) {
    padding-right: 0;
  }
}

@media (max-width: 575px) {
  .article-nav {
    flex-direction: column;
    align-items: stretch !important;
  }
  .article-nav-btn {
    justify-content: center;
    width: 100%;
  }
  .article-nav-empty {
    justify-content: center;
    width: 100%;
    margin: 0 !important;
  }
}
</style>
