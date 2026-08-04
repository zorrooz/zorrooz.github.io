<!-- Article.vue -->
<template>
  <div class="container view-container article-view">
    <div class="reading-progress" :style="{ width: progressPercent + '%' }" aria-hidden="true"></div>

    <div class="row py-4 px-0">
      <div class="col-12 col-lg-2 order-2 order-lg-1 docs-sidebar-col">
        <div class="sticky-sidebar" ref="leftSidebarContent">
          <div v-if="isDesktop" class="navigation-container mb-0">
            <NavigationTree />
          </div>
        </div>
      </div>

      <div class="col-12 col-lg-8 order-1 order-lg-2 docs-main-col" ref="mainContent">
        <div class="article-panel">
          <div class="article-panel__body">
            <div class="article-content">
              <div v-if="currentPost" class="article-meta">
                <h1 class="article-title">{{ currentPost.title }}</h1>
                <div class="article-meta__row">
                  <span v-if="isNote && currentPost.date" class="meta-line"><i class="fas fa-calendar-alt"></i>{{ updatedAtText }} {{ currentPost.date }}</span>
                  <span v-if="readingMinutes > 0" class="meta-line"><i class="fas fa-clock"></i>{{ getReadingTimeText(readingMinutes) }}</span>
                  <button type="button" class="article-copy-btn"
                    :class="{ 'article-copy-btn--copied': copyFeedback }" @click="copyArticle"
                    :aria-label="t('copyArticle')" aria-live="polite">
                    <i :class="copyFeedback ? 'fas fa-check' : 'fas fa-copy'"></i>
                    <span>{{ copyFeedback ? t('copied') : t('copyArticle') }}</span>
                  </button>
                </div>
                <div v-if="currentPost.tags?.length" class="article-meta__tags">
                  <span v-for="(tag, idx) in currentPost.tags" :key="idx" class="article-tag" @click="onTagClick(tag)">
                    # {{ tag }}
                  </span>
                </div>
              </div>

              <RenderMarkdown v-if="rawMarkdown" :rawMarkdown="rawMarkdown" :articlePath="currentPost?.path || ''"
                :articleTitle="currentPost?.title || ''" @markdown-rendered="handleMarkdownRendered" />

              <nav class="article-navigation" v-if="rawMarkdown">
                <span v-if="!prevPost && nextPost" class="article-nav-spacer" aria-hidden="true"></span>

                <router-link v-if="prevPost" :to="toArticle(prevPost.path)" class="article-nav-item prev">
                  <div class="nav-arrow"><i class="fas fa-arrow-left"></i></div>
                  <div class="nav-details">
                    <div class="nav-label">{{ prevPageText }}</div>
                    <div class="nav-title">{{ prevPost.title }}</div>
                  </div>
                </router-link>

                <router-link v-if="nextPost" :to="toArticle(nextPost.path)" class="article-nav-item next">
                  <div class="nav-details">
                    <div class="nav-label">{{ nextPageText }}</div>
                    <div class="nav-title">{{ nextPost.title }}</div>
                  </div>
                  <div class="nav-arrow"><i class="fas fa-arrow-right"></i></div>
                </router-link>
              </nav>
            </div>
          </div>
        </div>
      </div>

      <div class="col-12 col-lg-2 order-3 docs-toc-col">
        <div class="sticky-sidebar" ref="rightSidebarContent">
          <div v-if="isDesktop" class="toc-container mt-0">
            <OnThisPage ref="onThisPageRef" containerSelector=".markdown-body" :levels="[2, 3]" :offset="88" />
          </div>
        </div>
      </div>
    </div>

    <TocDrawer v-if="rawMarkdown" />
  </div>
</template>

<script setup lang="ts">
import { computed, nextTick, onBeforeUnmount, onMounted, ref, useTemplateRef, watch } from 'vue'
import { useI18n } from 'vue-i18n'
import { useRoute, useRouter } from 'vue-router'
import RenderMarkdown from '@/components/layout/RenderMarkdown.vue'
import OnThisPage from '@/components/layout/OnThisPage.vue'
import TocDrawer from '@/components/widgets/TocDrawer.vue'
import NavigationTree from '@/components/layout/NavigationTree.vue'
import { loadCategories, loadHtmlContent, loadMarkdownSource } from '@/utils/contentLoader'
import { toLocalePath } from '@/utils/localePath'
import { readingTimeMinutes } from '@/utils/readingTime'

defineOptions({
  name: 'ArticleView',
  head() {
    return {
      title: this.currentPost ? `${this.currentPost.title} - zorrooz’s blog` : 'zorrooz’s blog',
      meta: this.currentPost?.description
        ? [{ name: 'description', content: this.currentPost.description }]
        : [],
    }
  }
})

const { t, locale } = useI18n()
const route = useRoute()
const router = useRouter()

defineProps({ path: { type: [String, Array], default: '' } })

const rawMarkdown = ref('')
const currentPath = ref('')
const allArticles = ref<any[]>([])
const categoryList = ref<any[]>([])
const viewportWidth = ref((typeof window !== 'undefined' ? window.innerWidth : 1024))
const progressPercent = ref(0)
const copyFeedback = ref(false)
let copyFeedbackTimer: number | null = null
let scrollTicking = false

const onThisPageRef = useTemplateRef<InstanceType<typeof OnThisPage>>('onThisPageRef')
const leftSidebarContent = useTemplateRef<HTMLElement>('leftSidebarContent')
const rightSidebarContent = useTemplateRef<HTMLElement>('rightSidebarContent')

const isDesktop = computed(() => viewportWidth.value >= 992)
const isNote = computed(() => !!(currentPost.value?.path && currentPost.value.path.startsWith('notes/')))
const updatedAtText = computed(() => t('updatedAt'))
const prevPageText = computed(() => t('prevPage'))
const nextPageText = computed(() => t('nextPage'))

const currentPost = computed(() => {
  if (!currentPath.value) return null

  return allArticles.value.find(article => {
    const articlePath = article.path.replace(/\.md$/, '')
    const cp = currentPath.value.replace(/\.md$/, '')

    if (articlePath === cp) return true
    if (cp.endsWith('-en') && articlePath === cp.replace(/-en$/, '')) return true
    if (!cp.endsWith('-en') && articlePath === cp + '-en') return true

    return false
  })
})

const groupLinearArticles = computed(() => {
  if (!currentPost.value) return []
  const [type, group] = currentPost.value.path.replace(/\.md$/, '').split('/')
  const linear: { title: string; path: string }[] = []

  const pushFromUrl = (title: string, articleUrl: string) => {
    if (!articleUrl) return
    const parts = String(articleUrl).replace(/^\/+/, '').split('/')
    const i0 = parts[0] === 'article' ? 1 : 0
    const t = parts[i0], g = parts[i0 + 1]
    if (t !== type || g !== group) return
    const rest = parts.slice(i0 + 2)
    const pathNoExt = `${t}/${g}/${rest.join('/')}`
    linear.push({ title, path: `${pathNoExt}.md` })
  }

  if (Array.isArray(categoryList.value)) {
    for (const section of categoryList.value) {
      if (!Array.isArray(section.items)) continue
      for (const item of section.items) {
        if (item?.name !== group) continue
        if (Array.isArray(item.articles)) {
          item.articles.forEach((a: any) => pushFromUrl(a.title, a.articleUrl))
        }
        if (Array.isArray(item.categories)) {
          item.categories.forEach((cat: any) => {
            if (Array.isArray(cat.articles)) {
              cat.articles.forEach((a: any) => pushFromUrl(a.title, a.articleUrl))
            }
          })
        }
      }
    }
  }
  return linear
})

const currentLinearIndex = computed(() => {
  if (!currentPost.value) return -1
  return groupLinearArticles.value.findIndex(a => a.path.replace(/\.md$/, '') === currentPost.value.path.replace(/\.md$/, ''))
})

const prevPost = computed(() => {
  const idx = currentLinearIndex.value
  return idx > 0 ? groupLinearArticles.value[idx - 1] : null
})

const nextPost = computed(() => {
  const idx = currentLinearIndex.value
  const last = groupLinearArticles.value.length - 1
  return idx >= 0 && idx < last ? groupLinearArticles.value[idx + 1] : null
})

const readingMinutes = computed(() => {
  const wc = currentPost.value?.wordCount
  return typeof wc === 'number' ? readingTimeMinutes(wc) : 0
})

function getReadingTimeText(minutes: number) {
  return t('articleReadingTime', { minutes })
}

function toArticle(p: string) {
  return toLocalePath(`/article/${p.replace(/\.md$/, '')}`)
}

function onTagClick(tag: string) {
  if (!tag) return
  const prefix = locale.value === 'en-US' ? '/en' : '/zh'
  router.push({ path: `${prefix}/`, query: { tag, from: route.fullPath } })
}

async function copyArticle() {
  if (!currentPost.value) return
  // 优先复制完整 markdown 源（含 frontmatter，粘贴即得 .md 文件）；加载失败回退纯文本
  let text = ''
  try {
    const md = await loadMarkdownSource(currentPost.value.path)
    if (md) {
      text = md.trim()
    }
  } catch {
    // fall through to text extraction
  }
  if (!text) {
    const body = document.querySelector('.markdown-body')
    if (!body) return
    const clone = body.cloneNode(true) as HTMLElement
    clone
      .querySelectorAll('.heading-anchor, .code-block-header, .copy-button, .table-copy-btn')
      .forEach(el => el.remove())
    const bodyText = (clone.innerText || '').trim()
    text = `${currentPost.value.title}\n\n${bodyText}`
  }
  try {
    await navigator.clipboard.writeText(text)
  } catch (err) {
    console.error(t('copyFailed'), err)
    const textArea = document.createElement('textarea')
    textArea.value = text
    document.body.appendChild(textArea)
    textArea.select()
    document.execCommand('copy')
    document.body.removeChild(textArea)
  } finally {
    copyFeedback.value = true
    if (copyFeedbackTimer) clearTimeout(copyFeedbackTimer)
    copyFeedbackTimer = window.setTimeout(() => {
      copyFeedback.value = false
    }, 1200)
  }
}

function buildFromCategories() {
  const all: any[] = []

  const pushArticle = (artTitle: string, articleUrl: string, tags: any[] = [], dateStr: string = '', wordCount = 0) => {
    if (typeof articleUrl !== 'string' || !articleUrl.trim()) return
    const parts = articleUrl.replace(/^\/+/, '').split('/')
    const idxArticle = parts[0] === 'article' ? 1 : 0
    const t = parts[idxArticle], g = parts[idxArticle + 1]
    const restParts = parts.slice(idxArticle + 2)
    const subKey = restParts[0] || '__root__'
    const rest = restParts.join('/')
    const pathNoExt = `${t}/${g}/${rest}`

    const art = { title: artTitle, path: `${pathNoExt}.md`, date: dateStr || '', tags: Array.isArray(tags) ? tags : [], preview: '', category: `${t}/${g}/${subKey}`, wordCount: typeof wordCount === 'number' ? wordCount : 0 }
    all.push(art)
  }

  try {
    const categoryData = loadCategories()

    if (Array.isArray(categoryData)) {
      categoryData.forEach((section: any) => {
        if (!Array.isArray(section.items)) return
        section.items.forEach((item: any) => {
          const itemLatest = item?.stats?.latestDate || ''
          if (Array.isArray(item.articles)) {
            item.articles.forEach((a: any) => pushArticle(a.title, a.articleUrl, a?.tags || [], itemLatest, a?.wordCount))
          }
          if (Array.isArray(item.categories)) {
            item.categories.forEach((cat: any) => {
              const catLatest = cat?.stats?.latestDate || itemLatest || ''
              if (Array.isArray(cat.articles)) {
                cat.articles.forEach((a: any) => pushArticle(a.title, a.articleUrl, a?.tags || [], catLatest, a?.wordCount))
              }
            })
          }
        })
      })
    }

    allArticles.value = all
    categoryList.value = categoryData
  } catch {
    allArticles.value = []
    categoryList.value = []
  }
}

function onResize() {
  viewportWidth.value = window.innerWidth
  updateSidebarDimensions()
}

function onScroll() {
  if (scrollTicking) return
  scrollTicking = true
  requestAnimationFrame(() => {
    scrollTicking = false
    const doc = document.documentElement
    const max = doc.scrollHeight - window.innerHeight
    progressPercent.value = max > 0 ? Math.min(100, (window.scrollY / max) * 100) : 0
  })
}

function normalizeRoutePathParam(p: string | string[]) {
  return Array.isArray(p) ? p.join('/') : (typeof p === 'string' && p.length ? p : '')
}

function loadArticleContent() {
  rawMarkdown.value = ''
  try {
    const currentPathClean = normalizeRoutePathParam(route.params.path)

    const isEnglish = locale.value === 'en-US'

    let matchedPost = allArticles.value.find(article => {
      const articlePath = article.path.replace(/\.md$/, '')

      if (isEnglish) {
        return articlePath === currentPathClean ||
          articlePath === currentPathClean.replace(/-en$/, '') + '-en'
      } else {
        return articlePath === currentPathClean ||
          articlePath === currentPathClean.replace(/-en$/, '')
      }
    })

    if (!matchedPost) {
      matchedPost = allArticles.value.find(article => {
        const articlePath = article.path.replace(/\.md$/, '')
        const cleanCurrentPath = currentPathClean.replace(/-en$/, '')
        return articlePath === cleanCurrentPath ||
          articlePath === cleanCurrentPath + '-en'
      })
    }

    if (!matchedPost) throw new Error(`Article not found: ${currentPathClean}`)
    currentPath.value = matchedPost.path

    rawMarkdown.value = loadHtmlContent(currentPath.value)

    nextTick(() => {
      if (typeof window === 'undefined') return
      window.scrollTo({ top: 0, behavior: 'smooth' })
      updateSidebarDimensions()
      onThisPageRef.value?.refreshToc()
    })
  } catch {
    rawMarkdown.value = '# Article Not Found\n\nThe requested article could not be loaded. Please check the URL.'
    nextTick(() => {
      if (typeof window === 'undefined') return
      window.scrollTo({ top: 0, behavior: 'smooth' })
      onThisPageRef.value?.refreshToc()
    })
  }
}

function updateSidebarDimensions() {
  if (typeof window === 'undefined') return
  const header = document.querySelector('header')
  const headerH = header?.offsetHeight || 60
  const viewportH = window.innerHeight
  const availableH = Math.max(200, viewportH - headerH - 24 - 24)
  const sidebarEls = [leftSidebarContent.value, rightSidebarContent.value]
  sidebarEls.forEach(el => {
    if (!el) return
    el.style.maxHeight = `${availableH}px`
    el.style.overflowY = 'auto'
  })
}

function handleMarkdownRendered() {
  updateSidebarDimensions()
  nextTick(() => onThisPageRef.value?.refreshToc())
}

buildFromCategories()
loadArticleContent()

watch(locale, (newLocale, oldLocale) => {
  if (newLocale !== oldLocale) {
    buildFromCategories()
    loadArticleContent()
  }
})

watch(() => route.params.path, (newPathParam, oldPathParam) => {
  const oldPath = normalizeRoutePathParam(oldPathParam)
  const newPath = normalizeRoutePathParam(newPathParam)
  if (oldPath !== newPath) {
    onThisPageRef.value?.resetToc()
    loadArticleContent()
  }
})

watch(rawMarkdown, () => {
  nextTick(() => updateSidebarDimensions())
})

onMounted(() => {
  viewportWidth.value = window.innerWidth
  window.addEventListener('resize', onResize)
  window.addEventListener('scroll', onScroll, { passive: true })
})

onBeforeUnmount(() => {
  window.removeEventListener('resize', onResize)
  window.removeEventListener('scroll', onScroll)
  if (copyFeedbackTimer) clearTimeout(copyFeedbackTimer)
})
</script>

<style scoped>
.article-view {
  max-width: 1280px;
}

.article-panel {
  background: transparent;
  border: none;
  border-radius: 0;
  box-shadow: none;
}

.reading-progress {
  position: fixed;
  top: 0;
  left: 0;
  height: 2px;
  width: 0;
  background: var(--primary);
  z-index: 1045;
  border-radius: 0 99px 99px 0;
  pointer-events: none;
  transition: width 0.08s linear;
}

.article-panel__body {
  padding: var(--sp-12) 0;
}

.sticky-sidebar {
  position: sticky;
  top: 92px;
  box-sizing: border-box;
  width: 100%;
  -webkit-overflow-scrolling: touch;
  transition: max-height 0.2s ease;
}

.article-content {
  min-height: 400px;
  max-width: 820px;
  margin: 0 auto;
}

.article-meta {
  margin-bottom: var(--sp-10);
  padding-bottom: var(--sp-8);
  border-bottom: 1px solid var(--line);
}

.article-meta__row {
  display: flex;
  flex-wrap: wrap;
  gap: var(--sp-4);
  align-items: center;
  margin-top: var(--sp-5);
  font-size: var(--text-sm);
  color: var(--fg-3);
}

.article-meta .meta-line {
  display: inline-flex;
  align-items: center;
  gap: 7px;
}

.article-meta .meta-line i {
  font-size: 12px;
  color: var(--primary-muted);
}

.article-meta__tags {
  display: flex;
  flex-wrap: wrap;
  gap: var(--sp-2);
  margin-top: var(--sp-5);
}

.article-copy-btn {
  display: inline-flex;
  align-items: center;
  gap: 6px;
  margin-left: auto;
  font-size: var(--text-sm);
  font-weight: 500;
  color: var(--fg-3);
  background: transparent;
  border: 1px solid transparent;
  border-radius: var(--radius-btn);
  padding: 4px 10px;
  cursor: pointer;
  transition: color 0.14s ease, background-color 0.14s ease;
}

.article-copy-btn:hover {
  color: var(--primary);
  background: var(--tint);
}

.article-copy-btn--copied {
  color: var(--primary);
}

.article-tag {
  font-size: var(--text-sm);
  font-weight: 500;
  color: var(--fg-2);
  background: var(--surface-2);
  padding: 4px 13px;
  border-radius: var(--radius-pill);
  transition: color 0.14s ease, background-color 0.14s ease;
  cursor: pointer;
}

.article-tag:hover {
  color: var(--primary);
  background: var(--tint);
}

@media (min-width: 992px) {
  .docs-sidebar-col {
    border-right: 1px solid var(--line);
    padding-right: var(--sp-5);
  }

  .docs-toc-col {
    padding-left: var(--sp-5);
  }

  .docs-main-col {
    padding-left: var(--sp-8);
    padding-right: var(--sp-8);
  }
}

:deep(.markdown-body) {
  font-size: var(--text-body);
  line-height: 1.75;
  color: var(--fg);
}

:deep(.markdown-body p) {
  margin-bottom: 1.1rem;
}

:deep(.markdown-body h1) {
  font-size: var(--text-4xl);
  font-weight: 700;
  margin-top: 3rem;
  margin-bottom: 1.25rem;
  letter-spacing: -0.02em;
}

:deep(.markdown-body h2) {
  font-size: var(--text-3xl);
  font-weight: 700;
  margin-top: 2.75rem;
  margin-bottom: 1rem;
  letter-spacing: -0.015em;
  padding-bottom: 0.5rem;
  border-bottom: 1px solid var(--line);
}

:deep(.markdown-body h3) {
  font-size: var(--text-xl);
  font-weight: 600;
  margin-top: 2.25rem;
  margin-bottom: 0.875rem;
  letter-spacing: -0.01em;
}

:deep(.markdown-body h4) {
  font-size: 18px;
  font-weight: 600;
  margin-top: 1.75rem;
  margin-bottom: 0.75rem;
}

:deep(.markdown-body h5) {
  font-size: var(--text-body);
  font-weight: 600;
  margin-top: 1.5rem;
  margin-bottom: 0.625rem;
}

:deep(.markdown-body h6) {
  font-size: var(--text-sm);
  font-weight: 600;
  margin-top: 1.25rem;
  margin-bottom: 0.5rem;
  color: var(--fg-2);
}

.article-navigation {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: var(--sp-4);
  margin-top: var(--sp-16);
  padding-top: var(--sp-8);
  border-top: 1px solid var(--line);
}

.article-nav-spacer {
  display: block;
}

.article-nav-item {
  display: flex;
  align-items: center;
  gap: var(--sp-4);
  padding: var(--sp-5) var(--sp-6);
  text-decoration: none;
  color: var(--fg);
  transition: color 0.14s ease, background-color 0.14s ease, border-color 0.14s ease,
    box-shadow 0.18s ease;
  min-width: 0;
  border: 1px solid var(--line);
  border-radius: var(--radius);
  background: var(--surface);
}

.article-nav-item:hover {
  border-color: color-mix(in srgb, var(--primary) 35%, transparent);
  background: var(--surface-2);
  box-shadow: var(--shadow-soft);
}

.article-nav-item.next {
  justify-content: flex-end;
  text-align: left;
}

.nav-arrow {
  font-size: 14px;
  color: var(--fg-3);
  transition: color 0.14s ease;
  flex-shrink: 0;
  display: flex;
  align-items: center;
  justify-content: center;
}

.article-nav-item:hover .nav-arrow {
  color: var(--primary);
}

.article-nav-item.next:hover .nav-arrow {
  transform: none;
}

.nav-details {
  display: flex;
  flex-direction: column;
  overflow: hidden;
  min-width: 0;
  gap: 3px;
}

.article-nav-item.next .nav-details {
  align-items: flex-end;
  flex: 1 1 auto;
}

.nav-label {
  font-size: var(--text-xs);
  color: var(--fg-3);
  transition: color 0.14s ease;
}

.article-nav-item:hover .nav-label {
  color: var(--primary);
}

.nav-title {
  font-size: var(--text-md);
  font-weight: 600;
  white-space: nowrap;
  overflow: hidden;
  text-overflow: ellipsis;
  direction: ltr;
  max-width: 300px;
}

.article-nav-item.next .nav-title {
  text-align: left;
}

@media (max-width: 991px) {
  .sticky-sidebar {
    position: static;
    top: auto;
    bottom: auto !important;
    max-height: none !important;
    overflow-y: visible !important;
  }
  .navigation-container {
    margin-bottom: 1rem;
    margin-top: 1rem;
  }

  .toc-container {
    margin-top: 1rem;
    margin-bottom: 1rem;
  }

  .article-panel__body {
    padding: var(--sp-8) 0;
  }
}

@media (max-width: 768px) {
  .article-view .row {
    --bs-gutter-x: 2.5rem;
  }

  .article-navigation {
    grid-template-columns: 1fr;
    gap: var(--sp-3);
  }

  .article-nav-spacer {
    display: none;
  }

  .article-nav-item.next {
    justify-content: flex-start;
    text-align: left;
  }

  .article-nav-item.next .nav-details {
    align-items: flex-start;
  }

  .article-nav-item.next .nav-title {
    text-align: left;
  }

  :deep(.markdown-body h1) {
    font-size: 26px;
  }

  :deep(.markdown-body h2) {
    font-size: 22px;
  }

  :deep(.markdown-body h3) {
    font-size: 19px;
  }
}

@media (max-width: 576px) {
  :deep(.markdown-body) {
    font-size: 16px;
  }

  :deep(.markdown-body h1) {
    font-size: 24px;
  }

  :deep(.markdown-body h2) {
    font-size: 20px;
  }

  :deep(.markdown-body h3) {
    font-size: 18px;
  }
}
</style>
