import type { AppMessages } from '../schema'

export default {
  // 通用
  tags: 'Tags',
  articles: 'Articles',
  words: 'Words',
  prevPage: 'Previous',
  nextPage: 'Next',
  pagination: 'Pagination',
  search: 'Search',
  theme: 'Theme',
  language: 'Language',
  menu: 'Menu',
  close: 'Close',
  searchPlaceholder: 'Search articles by title or content...',
  searchNoResults: 'No matching articles found',
  searchUnavailable: 'Search index failed to load',
  searchEscHint: 'Esc to close',
  searchResultsLabel: 'results',

  // 导航
  categories: 'Categories',
  resources: 'Resources',
  about: 'About',

  // 首页
  greeting: "Hello, I'm",
  greetingPrefix: '//',
  developer: 'Developer',
  wordUnit: 'words',
  recentPosts: 'Recent Posts',
  noPosts: 'No posts yet',

  // 分类页面
  notes: 'Notes',
  projects: 'Projects',
  topics: 'Topics',
  seeMore: 'See More',
  countPosts: '{count} posts',
  countProjects: '{count} projects',
  countTopics: '{count} topics',
  countWords: '{count} words',
  copyCode: 'Copy code',
  copyFailed: 'Copy failed',
  copyArticle: 'Copy article',
  copyTable: 'Copy table',
  copied: 'Copied',
  anchorHeading: 'Anchor to heading',
  // 未分类
  uncategorized: 'Uncategorized',

  // 文章页面
  updatedAt: 'Updated at',
  postReadingTime: '{minutes} min',
  articleReadingTime: 'Reading about {minutes} minutes',
  tableOfContents: 'On this page',
  openToc: 'Open table of contents',
  backToTop: 'Back to top',

  // 资源页面
  resourceSubtitle: 'Commonly used tools in bioinformatics and structural biology',

  // 404
  pageNotFound: 'This page could not be found',
  backHome: 'Back to home',

  // SEO title
  metaHome: 'Home',
  metaCategories: 'Categories',
  metaResources: 'Resources',
  metaAbout: 'About',
  metaNotFound: '404',

  // Accessibility
  skipToContent: 'Skip to content',
  articleNav: 'Article Navigation',

  // 关于页面
  experience: 'Experience',
  thanks: 'Thank you for your attention!',
  designedByPrefix: 'Designed by',
  designedBySuffix: '',

  // 过滤功能
  clearFilter: 'Clear filter',
  backToArticle: 'Back to article',
} satisfies AppMessages
